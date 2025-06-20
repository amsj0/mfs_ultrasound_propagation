import numpy as np
import sys

from scipy.signal import (convolve,convolve2d)
from util.heuristic import heuristic
from util.h5py_util import *
from util.store import Store
#from threading import Thread
from multiprocessing import Process,Lock,JoinableQueue
from os import cpu_count,remove
from config import parse_config, create_configfile

DELETE_MFS = False

def worker(q,a,l):
    """The process will continually pull elements from the shared queue 
    to process until reaching a None sentinel.
    """

    while True:

        current = q.get(timeout=5)  
        if current is None:
            q.task_done()
            break
        elt,MH,MR,path = current
        print("Starting to process elt: {}".format(elt))
        fn_analyse(elt,a,l,MH,MR,path)
  
        print("Finished processing elt: {}".format(elt))
        q.task_done()

def fn_analyse(elt,apod,lock,MH,MR,path):
    
    output_path,dataroot,heurisset = path
    
    datafile = dataroot + '_' + str(elt+1) + heurisset
    dataset  = dataroot + heurisset

    apod2 = apod[:,np.newaxis] @ apod[np.newaxis,:]

    domain = convolve(MH,apod[np.newaxis,:],mode='valid')
    response = convolve2d(MR,apod2,mode='valid')        
    print('DataFile {} read'.format(datafile))

    
    lock.acquire()
    # Do not use a proxy object from more than one thread unless you protect it with a lock.
    try:
        append_keyvalue_to_hdf5('doma', domain, elt, output_path + 'doma_', dataset)
        append_keyvalue_to_hdf5('resp', response, elt, output_path + 'resp_', dataset)
    finally:
        lock.release()

    print('DataFile {} appended'.format(output_path + '_' + dataset))


def analyse(name,store,config_file, output_path):

    config_tuple = create_configfile(parse_config,config_file, output_path)

    T,S,D,R,Neltoverlambda,nRD,g = config_tuple

    ppt_per_surface = 1 + int(np.floor(g.piston_radius * (Neltoverlambda / 100)))

    k0, kr, kr_length, dr, dr_length, sfr, RD, lambda0, d_cur = heuristic(nRD, g)

    dataroot = g.convergemod + '_' + str(g.nff) + '_' + str(int(g.iff*g.model_scale*100)) + '_' + str(int(g.fff*g.model_scale*100)) 

    rshape = (R.c.size-ppt_per_surface+1,T.c.size-ppt_per_surface+1)
    dshape = (D.c.size,T.c.size-ppt_per_surface+1)

    domaset_size = (sfr.size,) + dshape
    respset_size = (sfr.size,) + rshape

    # Ensure range_ppt has ppt_per_surface points, centered around 0
    half = ppt_per_surface // 2
    if ppt_per_surface % 2 == 0:
        range_ppt = np.arange(-half, half) + 0.5
    else:
        range_ppt = np.arange(-half, half + 1)

    apod = 8/(ppt_per_surface*np.pi)*np.sqrt(int(ppt_per_surface/2)**2-range_ppt**2)
    # apod = np.ones((ppt_per_surface,))
    
    processes = []

    cpc = cpu_count()

    lock = Lock()
    
    path = [output_path,dataroot,'_']

    queue_max_size = 20
    q = JoinableQueue(maxsize=queue_max_size)
    # Create processes
    processes = []
    num_proc = min(g.ffu-g.ifu,cpc)
    
    for i in range(num_proc):
        p = Process(target=worker, args=(q,apod,lock))
        p.start()
        processes.append(p)

    for jj in range(kr_length):
        for pp in range(dr_length):
            path[2] = '_' + str(int(g.skr[jj])) + '_' + str(int(g.sdr[pp]))

            dataset = path[1] + path[2]

            create_keysized_to_hdf5('doma', domaset_size, output_path + 'doma_', dataset)
            create_keysized_to_hdf5('resp', respset_size, output_path + 'resp_', dataset)

            for elt in range(g.ifu-1,g.ffu):
                
                datafile = path[1] + '_' + str(elt+1) + path[2]

                if '' in store.Solution:
                    with h5py.File(path[0] + datafile + '.h5', 'r') as f:
                        MH = np.copy(f['domain'])
                        MR = np.copy(f['receiver'])

                    if DELETE_MFS:
                        remove(path[0] + datafile + '.h5')
                else:
                    MH = store.Solution[datafile]['domain']
                    MR = store.Solution[datafile]['receiver']

                q.put((elt,MH,MR,path), block=True)

    # Add sentinels to shut down processes.
    print('Adding sentinels...')
    for _ in range(num_proc):
        q.put(None)

    q.join() # Block until all elements have been processed.

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        raise ValueError('Invalid number of arguments. Usage: {} input_file.yaml /path/to/output'.format(sys.argv[0]))

    input_file = sys.argv[1]
    output_path = sys.argv[2]

    store = Store('')
    analyse("analyse",store,input_file, output_path)