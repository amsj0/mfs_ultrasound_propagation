"""
Concurrent analysis of mesh-free solutions with apodization.
- Loads geometry from YAML config or Store.
- Discretizes geometry.
- Applies apodization via convolution.
- Saves results to HDF5 or Store object.

author: Agesinaldo Silva
date: June, 2024
"""

from __future__ import annotations

import sys
from concurrent.futures import ProcessPoolExecutor, as_completed
from os import cpu_count, remove
from typing import Iterable, Tuple
import numpy as np
import h5py

from util.heuristic import heuristic
from util.h5py_util import create_keysized_to_hdf5, append_keyvalue_to_hdf5
from util.store import Store
from config import parse_config, create_configfile

# Optional: use build_geometry wrapper if present
try:
    from config import build_geometry  # type: ignore
except Exception:
    def build_geometry(filename: str, output_path: str):
        return create_configfile(parse_config, filename, output_path)

DELETE_MFS = False


# ------------------------------- Worker --------------------------------------

def compute_analyse(elt: int, apod: np.ndarray, MH: np.ndarray, MR: np.ndarray):
    """
    Pure function executed in a worker process.
    Parameters
    ----------
    elt : int
        Frequency/element index.
    apod : np.ndarray, shape (P,)
        Apodization vector over piston-taps-per-surface (PPT).
    MH : np.ndarray
        Domain (mesh-free) slice for this elt.
    MR : np.ndarray
        Receiver slice for this elt.

    Returns
    -------
    domain : np.ndarray
        Convolution of MH with apod along transmitter axis (valid mode).
    respon : np.ndarray
        2D convolution of MR with outer(apod, apod) (valid mode).
    """
    apod2  = apod[:, np.newaxis] @ apod[np.newaxis, :]
    from scipy.signal import convolve, convolve2d
    domain = convolve(MH, apod[np.newaxis, :], mode='valid')
    respon = convolve2d(MR, apod2, mode='valid')
    return domain, respon


# ------------------------------ Helpers --------------------------------------

def load_config_or_store(store: Store, input_file: str, output_path: str):
    """
    Decide whether to load geometry/config from disk (HDF5) or from the Store.
    Returns (T, S, D, R, Neltoverlambda, nRD, g, dataroot).
    """
    cfg_tuple = build_geometry(input_file, output_path)
    if '' in store.Solution:
        # Running from HDF5
        T, S, D, R, Neltoverlambda, nRD, g = cfg_tuple
        dataroot = f"{g.convergemod}_{g.nff}_{int(g.iff*2*g.model_scale*100)}_{int(g.fff*2*g.model_scale*100)}"
    else:
        # Running from in-memory store
        config = store.Configuration
        dataroot = list(config.keys())[0]
        g = config[dataroot]
        T, S, D, R, Neltoverlambda, nRD, _ = cfg_tuple
    return T, S, D, R, Neltoverlambda, nRD, g, dataroot


def build_apodization(ppt_per_surface: int) -> np.ndarray:
    """
    Construct the apodization vector used for piston-taps-per-surface (PPT).

    This preserves the legacy formula but explains even/odd handling:
    - Even PPT: symmetric half-steps around 0 (e.g., -1.5, -0.5, 0.5, 1.5)
    - Odd  PPT: integer steps including 0 (e.g., -2, -1, 0, 1, 2)
    """
    half = ppt_per_surface // 2
    if ppt_per_surface % 2 == 0:
        rng = np.arange(-half, half) + 0.5
    else:
        rng = np.arange(-half, half + 1)
    return 8/(ppt_per_surface*np.pi) * np.sqrt(int(ppt_per_surface/2)**2 - rng**2)


def plan_jobs(g, kr_length: int, dr_length: int):
    """Yield all (jj, pp, elt) job triplets in sweep order."""
    for jj in range(kr_length):
        for pp in range(dr_length):
            for elt in range(g.ifu - 1, g.ffu):
                yield jj, pp, elt


def read_mfs_slice(store: Store, output_path: str, datafile: str):
    """Read MH/MR arrays for a given datafile, from HDF5 or the Store."""
    if '' in store.Solution:
        with h5py.File(output_path + datafile + '.h5', 'r') as f:
            MH = np.copy(f['domain'])
            MR = np.copy(f['receiver'])
        if DELETE_MFS:
            remove(output_path + datafile + '.h5')
    else:
        MH = store.Solution[datafile]['domain']
        MR = store.Solution[datafile]['receiver']
    return MH, MR


# ------------------------------ Main API -------------------------------------

def analyse_concurrent(store: Store, input_file: str, output_path: str, task_name: str) -> Store:
    """
    Convolve mesh-free solutions with the apodization window concurrently and
    append results into HDF5 (or the in-memory Store).

    Parameters
    ----------
    store : Store
        In-memory storage; if empty ('' key present), reads inputs from HDF5.
    input_file : str
        YAML config path.
    output_path : str
        Directory to read .h5 inputs and write results.
    task_name : str
        'analyse' or 'mfsolution_analyse' -> write to HDF5;
        else -> write back into Store.
    """
    # 1) Load configuration and geometry
    T, S, D, R, Neltoverlambda, nRD, g, dataroot = load_config_or_store(store, input_file, output_path)

    # 2) Heuristics and shapes
    k0, kr, kr_length, dr, dr_length, sfr, RD, lambda0, d_cur = heuristic(nRD, g)
    ppt_per_surface = 1 + int(np.floor(g.piston_radius * (Neltoverlambda / 100)))

    rshape = (R.c.size - ppt_per_surface + 1, T.c.size - ppt_per_surface + 1)
    dshape = (D.c.size,                      T.c.size - ppt_per_surface + 1)
    domaset_size  = (sfr.size,) + dshape
    respset_size  = (sfr.size,) + rshape

    # 3) Apodization vector
    apod = build_apodization(ppt_per_surface)

    # 4) Prepare outputs per (jj, pp)
    for jj in range(kr_length):
        for pp in range(dr_length):
            heurisset = f"_{int(g.skr[jj])}_{int(g.sdr[pp])}"
            dataset   = dataroot + heurisset
            if task_name in ("analyse", "mfsolution_analyse"):
                create_keysized_to_hdf5('doma', domaset_size, output_path + 'doma_', dataset)
                create_keysized_to_hdf5('resp', respset_size, output_path + 'resp_', dataset)
            else:
                store.init_store_behavior('doma_' + dataset, domaset_size)
                store.init_store_behavior('resp_' + dataset, respset_size)

    # 5) Build process pool
    max_procs = max(1, (cpu_count() or 1))
    procs = min(max_procs, max(1, g.ffu - g.ifu + 1))

    submitted = 0
    futures = {}

    with ProcessPoolExecutor(max_workers=procs) as ex:
        # Submit all jobs
        for jj, pp, elt in plan_jobs(g, kr_length, dr_length):
            heurisset = f"_{int(g.skr[jj])}_{int(g.sdr[pp])}"
            dataset   = dataroot + heurisset
            datafile  = f"{dataroot}_{elt+1}{heurisset}"

            MH, MR = read_mfs_slice(store, output_path, datafile)
            fut = ex.submit(compute_analyse, elt, apod, MH, MR)
            futures[fut] = (elt, dataset)
            submitted += 1

        # Collect results and write
        completed = 0
        for fut in as_completed(futures):
            elt, dataset = futures.pop(fut)
            domain, respon = fut.result()

            if task_name in ("analyse", "mfsolution_analyse"):
                append_keyvalue_to_hdf5('doma', domain, elt, output_path + 'doma_', dataset)
                append_keyvalue_to_hdf5('resp', respon, elt, output_path + 'resp_', dataset)
            else:
                result = {'doma_' + dataset: domain, 'resp_' + dataset: respon}
                store.load_dict_to_store_behaviour(result, elt)

            if (completed := completed + 1) % 50 == 0 or completed == submitted:
                print(f"Processed {completed}/{submitted} tasks")

    return store


# --------------------------------- CLI ---------------------------------------

if __name__ == "__main__":
    if len(sys.argv) != 3:
        raise ValueError(f"Invalid number of arguments. Usage: {sys.argv[0]} input_file.yaml /path/to/output")
    input_file = sys.argv[1]
    output_path = sys.argv[2]

    store = Store('','')
    analyse_concurrent(
        store=store,
        input_file=input_file,
        output_path=output_path,
        task_name="analyse"
    )