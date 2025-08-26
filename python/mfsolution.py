import numpy as np
import sys
from concurrent.futures import ThreadPoolExecutor, as_completed
import copy
from time import sleep
from util.store import Store
from compute import Compute
from config import *
from util.heuristic import heuristic
from util.h5py_util import *

class structtype():
    pass

def reconfigure(config_tuple):

    T_old,S_old,D_old,R_old,_,nRD,g = config_tuple

    S = {
        'lower': structtype(),
        'upper': structtype(),
        'collo': structtype()
    }
        
    S['lower'].c = S_old.ci
    S['upper'].c = S_old.co

    S['collo'].c = S_old.c
    S['collo'].n = S_old.n
    S['collo'].a = S_old.a

    T = structtype()

    T = {
        'emitter': structtype()
    } 
    
    T['emitter'].c = T_old.c
    T['emitter'].n = T_old.n
    T['emitter'].a = T_old.a
    T['emitter'].m = T_old.ndx

    P = {
        'domain': structtype(),
        'receiver': structtype()
    }

    P['domain'].c = D_old.c
    P['domain'].m = D_old.ndx

    P['receiver'].c = R_old.c
    P['receiver'].m = R_old.ndx

    return T,P,S,nRD,g

# Helper to run both lower-side and upper-side in the background
def _compute_full(P, T, S, p_out, p_cur):
    cp = Compute(P, T, S)
    cp.InitCL("GPU")
    cp.compute_lower_side(p_out)
    cp.compute_upper_side(p_cur, p_out)
    return cp

# Worker to run scatter propagation and saving in background
# Delays propagate_scatter until inside save thread
def _scatter_and_save(cp_inst, datafile, outpath, task_name, ST):
    M = cp_inst.propagate_scatter()
    if task_name == "__main__":
        ST.load_dict_to_store(M, datafile)
    else:
        save_dict_to_hdf5(M, outpath, datafile)

def mfsolution(name, config_file, output_path):
    print(f"Config: {config_file}")
    print(f"Output path: {output_path}")
    T, P, S, nRD, g = reconfigure(create_configfile(parse_config, config_file, output_path))

    k0, kr, kr_length, dr, dr_length, sfr, RD, lambda0, d_cur = heuristic(nRD, g)
    dataroot = f"{g.convergemod}_{g.nff}_{int(g.iff*g.model_scale*100)}_{int(g.fff*g.model_scale*100)}"

    ST = Store('')

    # Pools: one for compute, one for saving
    compute_pool = ThreadPoolExecutor(max_workers=1)
    save_pool = ThreadPoolExecutor(max_workers=1)
    
    # Aggregate all compute futures across ii, jj, pp
    compute_futures = {}
    try:
        for ii in range(g.ifu - 1, g.ffu):
            # Prepare parameters
            k_cur = k0[ii] * lambda0 / RD
            p_out = [k_cur, d_cur]

            for jj in range(kr_length):
                for pp in range(dr_length):
                    dc_cur = g.rjR(kr[jj], dr[pp])[ii]
                    kc_cur = g.keq(kr[jj], dr[pp])[ii]
                    p_cur = [kc_cur, dc_cur]
                    datafile = (
                        f"{dataroot}_{ii+1}_"
                        f"{int(g.skr[jj])}_{int(g.sdr[pp])}"
                    )
                    future = compute_pool.submit(
                        _compute_full, P, T, S, p_out, p_cur
                    )
                    compute_futures[future] = (datafile, name, output_path)

        # As each compute future completes, handle propagate and schedule scatter+save
        print("Waiting for computations to complete...")
        for comp_future in as_completed(compute_futures):
                cp_inst = comp_future.result()
                datafile, task_name, outpath = compute_futures[comp_future]

                # Sequential propagate_transfer (heavy CPU)
                cp_inst.propagate_transfer()

                # Background scatter + save: propagate_scatter runs here
                save_future = save_pool.submit(
                    _scatter_and_save,
                    cp_inst,
                    datafile,
                    outpath,
                    task_name,
                    ST,
                )
                
                save_future.add_done_callback(
                    lambda f, df=datafile: print(f"Datafile {df} saved")
                )                
    finally:
        # Ensure pools shut down
        compute_pool.shutdown(wait=True)
        save_pool.shutdown(wait=True)

    return ST

if __name__ == "__main__":

    if len(sys.argv) != 3:
        raise ValueError(f'Invalid number of arguments. Usage: {sys.argv[0]} input_file.yaml')
    input_file, output_path = sys.argv[1], sys.argv[2]
    mfsolution("mfsolution", input_file, output_path)
