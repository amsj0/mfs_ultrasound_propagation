import numpy as np
import sys

from scipy.signal import (convolve,convolve2d)
from util.heuristic import heuristic
from util.h5py_util import *
from util.store import Store
#from threading import Thread
from concurrent.futures import ProcessPoolExecutor, as_completed
from os import cpu_count,remove
from config import parse_config, create_configfile

DELETE_MFS = False

def compute_analyse(elt, apod, MH, MR, path_tuple, task_name, dataset):
    """
    Pure function executed in a worker process.
    Returns (elt, (domain, respon), dataset).
    """
    output_path, dataroot, heurisset = path_tuple
    # Local (workers don't touch HDF5; parent will)
    apod2  = apod[:, np.newaxis] @ apod[np.newaxis, :]
    domain = convolve(MH, apod[np.newaxis, :], mode='valid')
    respon = convolve2d(MR, apod2, mode='valid')
    return elt, (domain, respon), dataset


def analyse_concurrent(
        store: Store, input_file: str, output_path: str, task_name: str) -> Store:


    config_tuple = create_configfile(parse_config, input_file, output_path)

    if '' in store.Solution:
        T, M, S, D, R, Neltoverlambda, nRD, g = config_tuple
        dataroot = f"{g.convergemod}_{g.nff}_{int(g.iff*g.model_scale*100)}_{int(g.fff*g.model_scale*100)}"
    else:
        config = store.Configuration
        dataroot = list(config.keys())[0]
        g = config[dataroot]
        T, M, S, D, R, Neltoverlambda, nRD, _ = config_tuple

    ppt_per_surface = 1 + int(np.floor(g.piston_radius * (Neltoverlambda / 100)))
    k0, kr, kr_length, dr, dr_length, sfr, RD, lambda0, d_cur = heuristic(nRD, g)

    rshape = (R.c.size - ppt_per_surface + 1, T.c.size - ppt_per_surface + 1)
    dshape = (D.c.size,                      T.c.size - ppt_per_surface + 1)

    domaset_size  = (sfr.size,) + dshape
    respset_size  = (sfr.size,) + rshape

    # Apodization
    half = ppt_per_surface // 2
    if ppt_per_surface % 2 == 0:
        range_ppt = np.arange(-half, half) + 0.5
    else:
        range_ppt = np.arange(-half, half + 1)
    apod = 8/(ppt_per_surface*np.pi) * np.sqrt(int(ppt_per_surface/2)**2 - range_ppt**2)

    # Build executor
    num_proc = min(g.ffu - g.ifu, (cpu_count() or 1))
    if num_proc <= 0:
        num_proc = 1

    # For reading per-task source data
    resp = None  # (not used here; we read MH/MR per file/task below)

    total_tasks = kr_length * dr_length * (g.ffu - g.ifu + 1)
    submitted = 0

    with ProcessPoolExecutor(max_workers=num_proc) as ex:
        futures = []

        # Submit all tasks
        for jj in range(kr_length):
            for pp in range(dr_length):
                heurisset = f"_{int(g.skr[jj])}_{int(g.sdr[pp])}"
                dataset   = dataroot + heurisset
                path      = [output_path, dataroot, heurisset]  # make a fresh (copy-like) object
                path_tuple = (path[0], path[1], path[2])

                # Create destination datasets once per (jj,pp)
                if task_name in ("analyse", "mfsolution_analyse"):
                    create_keysized_to_hdf5('doma', domaset_size, output_path + 'doma_', dataset)
                    create_keysized_to_hdf5('resp', respset_size, output_path + 'resp_', dataset)
                else:
                    store.init_store_behavior('doma_' + dataset, domaset_size)
                    store.init_store_behavior('resp_' + dataset, respset_size)

                for elt in range(g.ifu - 1, g.ffu):
                    datafile = f"{dataroot}_{elt+1}{heurisset}"

                    # Load MH/MR (parent reads; workers get arrays by pickle)
                    if '' in store.Solution:
                        with h5py.File(output_path + datafile + '.h5', 'r') as f:
                            MH = np.copy(f['domain'])
                            MR = np.copy(f['receiver'])
                        if DELETE_MFS:
                            remove(output_path + datafile + '.h5')
                    else:
                        MH = store.Solution[datafile]['domain']
                        MR = store.Solution[datafile]['receiver']

                    fut = ex.submit(compute_analyse, elt, apod, MH, MR, path_tuple, task_name, dataset)
                    futures.append(fut)
                    submitted += 1

        # Collect results as they complete; write outputs in the parent
        completed = 0
        for fut in as_completed(futures):
            elt, (domain, respon), dataset = fut.result()

            if task_name in ("analyse", "mfsolution_analyse"):
                # Parent writes to HDF5 (no locks needed)
                append_keyvalue_to_hdf5('doma', domain, elt, output_path + 'doma_', dataset)
                append_keyvalue_to_hdf5('resp', respon, elt, output_path + 'resp_', dataset)
            else:
                result = {'doma_' + dataset: domain, 'resp_' + dataset: respon}
                store.load_dict_to_store_behaviour(result, elt)

            completed += 1
            if completed % 50 == 0 or completed == submitted:
                print(f"Processed {completed}/{submitted} tasks")

    return store


if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        raise ValueError('Invalid number of arguments. Usage: {} input_file.yaml /path/to/output'.format(sys.argv[0]))

    input_file = sys.argv[1]
    output_path = sys.argv[2]

    store = Store('','')
    
    analyse_concurrent(
        store=store,
        input_file=input_file,
        output_path=output_path,
        task_name="analyse"
    )