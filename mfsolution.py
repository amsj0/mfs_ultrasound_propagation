""""
Mesh-free solution pipeline.
- Loads geometry from YAML config.
- Discretizes geometry.
- Runs mesh-free solution across frequency and ratio sweeps.
- Saves results to HDF5 or Store object.

author: Agesinaldo Silva
date: June, 2024
"""

from __future__ import annotations

import sys
import logging
import concurrent.futures as cf
#from concurrent.futures import ThreadPoolExecutor, as_completed
from typing import Tuple, Dict, Any, Callable, Iterable
import numpy as np
import gc, weakref

from util.store import Store
from compute import Compute
from util.heuristic import heuristic
from util.h5py_util import save_dict_to_hdf5, load_hdf5_to_array, save_dict_to_hdf5_toeplitz

# Import config helpers explicitly (no star import)
from config import parse_config, create_configfile

# ------------------------------- Setup ---------------------------------------

class structtype:
    """Simple mutable container, compatible with the rest of the codebase."""
    pass


def _cleanup_cp(cp: Compute) -> None:
    """Best-effort cleanup of Compute / GPU resources."""
    try:
        for meth in ("release",):
            if hasattr(cp, meth):
                getattr(cp, meth)()
    finally:
        del cp
        gc.collect()


# ----------------------------- Configuration ---------------------------------

def load_geometry(config_file: str, output_path: str):
    """
    Load and discretize the full geometry given a YAML configuration.

    Returns
    -------
    T, M, S, D, R, Neltoverlambda, nRD, g
        Objects from `create_configfile`, and the parsed global parameters `g`.

    Notes
    -----
    - Internally calls `parse_config(config_file)` to build `g`.
    - Then calls `create_configfile(parse_config, ...)` which both returns
      the discretized geometry and saves it to HDF5 for reproducibility.
    """
    return create_configfile(parse_config, config_file, output_path)


def map_geometry_to_compute_inputs(
    T_src, S_src, D_src, R_src
) -> Tuple[Dict[str, Any], Dict[str, Any], Dict[str, Any]]:
    """
    Map raw geometry objects (T/S/D/R) to the dictionary structures expected
    by `Compute(P, T, S)`.

    Returns
    -------
    T_map : dict
        {'emitter': ..., 'mirror': ...}
    P_map : dict
        {'domain': ..., 'receiver': ...}
    S_map : dict
        {'lower': ..., 'upper': ..., 'collo': ...}
    """
    # Surrounding surface (S): three views used by the solver
    S_map = {
        'lower': structtype(),
        'upper': structtype(),
        'collo': structtype(),
    }
    S_map['lower'].c = getattr(S_src, 'ci')
    S_map['upper'].c = getattr(S_src, 'co')
    S_map['collo'].c = getattr(S_src, 'c')
    S_map['collo'].n = getattr(S_src, 'n')
    S_map['collo'].a = getattr(S_src, 'a')

    # Transmitter/emitter (T) and mirror (M)
    T_map = {'emitter': structtype(), 'mirror': structtype()}
    T_map['emitter'].c = getattr(T_src, 'c')
    T_map['emitter'].n = getattr(T_src, 'n')
    T_map['emitter'].a = getattr(T_src, 'a')
    T_map['emitter'].m = getattr(T_src, 'ndx')

    # Problem/domain (D) and receiver (R)
    P_map = {'domain': structtype(), 'receiver': structtype()}
    P_map['domain'].c = getattr(D_src, 'c')
    P_map['domain'].m = getattr(D_src, 'ndx')

    P_map['receiver'].c = getattr(R_src, 'c')
    P_map['receiver'].m = getattr(R_src, 'ndx')

    return T_map, P_map, S_map


# ------------------------------ Workers --------------------------------------

def _compute_full(P_map, T_map, S_map, p_out, p_cur) -> Compute:
    cp = Compute(P_map, T_map, S_map)
    cp.InitCL("CPU", prefer_fp64=True)
    cp.compute_lower_side(p_out)
    cp.compute_upper_side(p_cur, p_out)
    return cp

def _load_propagate_matrix_from_datafile(cp_inst: Compute, datafile: str, outpath: str) -> None:
    MT = load_hdf5_to_array(datafile, outpath, 'MT')
    cp_inst.set_transfer(MT)

def _save_propagate_matrix_to_datafile(cp_inst: Compute, datafile: dict, outpath: str) -> None:
    save_dict_to_hdf5_toeplitz({'MT': cp_inst.MT}, outpath, datafile)

def _scatter_and_save(cp_inst: Compute, datafile: str, g, outpath: str, task_name: str, ST: Store) -> None:
    M = cp_inst.propagate_scatter()
    if task_name == "mfsolution":
        save_dict_to_hdf5(M, outpath, datafile)
    else:
        ST.load_dict_to_store_solution(M, datafile)


# ----------------------------- Main pipeline ---------------------------------
def bounded_submit(
    jobs: Iterable[tuple],
    worker: Callable[[tuple], Any],
    max_in_flight: int = 1
):
    """
    Keep at most `max_in_flight` futures running. Yields (job_tuple, result, exc)
    as each future completes. With max_in_flight=1, you guarantee only one
    _compute_full (and thus one cl.Context) at a time.
    """
    it = iter(jobs)
    with cf.ThreadPoolExecutor(max_workers=max_in_flight, thread_name_prefix="gpu") as ex:
        pending: Dict[cf.Future, tuple] = {}

        # prime one job
        first = next(it, None)
        if first is not None:
            pending[ex.submit(worker, *first[:-2])] = first

        while pending:
            done, _ = cf.wait(pending, return_when=cf.FIRST_COMPLETED)
            for fut in done:
                job_done = pending.pop(fut)
                res = fut.result()

                # *** PREFETCH NEXT BEFORE YIELD ***
                nxt = next(it, None)
                if nxt is not None:
                    pending[ex.submit(worker, *nxt[:-2])] = nxt

                yield {job_done[-1]: res,
                       "payload": job_done[-2],}


def mfsolution(config_file: str, output_path: str, task_name: str) -> Store:
    """
    Run the mesh-free solution across all (frequency, ratio) combinations
    specified in the YAML config, and save results per sweep.

    Parameters
    ----------
    task_name : str
        When "mfsolution", save to HDF5; otherwise load into the `Store` object.
    config_file : str
        Path to YAML configuration file.
    output_path : str
        Output directory.

    Returns
    -------
    Store
        A `Store` initialized with metadata and used for non-HDF5 saving mode.
    """
    print(f"Config: {config_file}")
    print(f"Output path: {output_path}")

    # 1) Load parsed parameters and discretized geometry
    T_src, S_src, D_src, R_src, _, nRD, g = load_geometry(config_file, output_path)

    # 2) Convert raw geometry into the structures expected by Compute
    T_map, P_map, S_map = map_geometry_to_compute_inputs(T_src, S_src, D_src, R_src)

    # 3) Precompute heuristics and naming
    k0, kr, kr_length, dr, dr_length, sfr, RD, lambda0, d_cur = heuristic(nRD, g)
    
    dataroot = f"{g.nff}_{int(g.iff*2*g.model_scale*100)}_{int(g.fff*2*g.model_scale*100)}"
    datamod = f"{g.convergemod}_{dataroot}"
    ST = Store('', datamod, g)

    # 4) Thread pools
    #compute_executor = cf.ThreadPoolExecutor(max_workers=1)
    save_executor = cf.ThreadPoolExecutor(max_workers=1, thread_name_prefix="saver")

    jobs = []
    try:
        for ii in range(g.ifu - 1, g.ffu):
            k_cur = k0[ii] * lambda0 / RD
            p_out = [k_cur, d_cur]

            for jj in range(kr_length):
                for pp in range(dr_length):
                    dc_cur = g.rjR(kr[jj], dr[pp])[ii]
                    kc_cur = g.keq(kr[jj], dr[pp])[ii]
                    p_cur = [kc_cur, dc_cur]
                    map = f"{ii+1}_{int(g.skr[jj])}_{int(g.sdr[pp])}"
                    root_map = f"{dataroot}_{map}"
                    mod_map = f"{datamod}_{map}"
                    jobs.append((P_map, T_map, S_map, p_out, p_cur,root_map,mod_map))

        print("Waiting for computations to complete...")
        for item in bounded_submit(jobs, _compute_full, max_in_flight=1):
            

            [(job_id, compute),job_payload] = item.items()

            # Attempt to load propagate matrix
            try:
                _load_propagate_matrix_from_datafile(compute, job_payload[-1], output_path)
            except Exception as e:
                logging.warning(f"Failed to load propagate matrix for {job_payload}: {e}")

            # Transfer
            if compute.MT_is_not_set:
                print(f"Transfer matrix loaded for {job_payload}, proceeding to transfer.")
                compute.propagate()
                _save_propagate_matrix_to_datafile(compute, job_payload[-1], output_path)

            compute.transfer()

            # Scatter + save
            save_future = save_executor.submit(
                _scatter_and_save, compute, job_id, g, output_path, task_name, ST
            )
            save_future.add_done_callback(lambda fut, df=job_id: print(f"Datafile {df} saved"))

            # Cleanup after save completes
            cp_ref = weakref.ref(compute)
            def _after_save(_fut, _r=cp_ref):
                cp = _r()
                if cp is not None:
                    _cleanup_cp(cp)
            save_future.add_done_callback(_after_save)

            # Drop strong ref here (the save task keeps it alive until done)
            compute = None
            job_id = None

    finally:
        #compute_executor.shutdown(wait=True)
        save_executor.shutdown(wait=True)

    return ST


# ------------------------------ CLI entry ------------------------------------

if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise ValueError(f"Invalid number of arguments. Usage: {sys.argv[0]} config.yaml /path/to/output")
    input_file, output_path = sys.argv[1], sys.argv[2]
    mfsolution(input_file, output_path, "mfsolution")