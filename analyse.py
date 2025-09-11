import numpy as np
import sys

from scipy.signal import (convolve,convolve2d)
from util.heuristic import heuristic
from util.h5py_util import *
from util.store import Store
#from threading import Thread
from multiprocessing import Process,Lock,JoinableQueue, Queue
from os import cpu_count,remove
from config import parse_config, create_configfile

DELETE_MFS = False

def worker(q,r,a,l):
    """The process will continually pull elements from the shared queue 
    to process until reaching a None sentinel.
    """

    while True:
        
        current = q.get(timeout=25)  
        if current is None:
            q.task_done()
            break
        elt,MH,MR,path,task_name,dataset = current
        print("Starting to process elt: {}".format(elt))
        result = fn_analyse(elt,a,l,MH,MR,path,task_name)
        r.put((elt,result,dataset))
        print("Finished processing elt: {}".format(elt))
        q.task_done()

def fn_analyse(elt,apod,lock,MH,MR,path,task_name):
    
    output_path,dataroot,heurisset = path
    
    datafile = dataroot + '_' + str(elt+1) + heurisset
    dataset  = dataroot + heurisset

    apod2 = apod[:,np.newaxis] @ apod[np.newaxis,:]

    domain = convolve(MH,apod[np.newaxis,:],mode='valid')
    respon = convolve2d(MR,apod2,mode='valid')        
    print('DataFile {} read'.format(datafile))

    if task_name == "analyse" or task_name == "mfsolution_analyse":
        lock.acquire()
        # Do not use a proxy object from more than one thread unless you protect it with a lock.
        try:
            append_keyvalue_to_hdf5('doma', domain, elt, output_path + 'doma_', dataset)
            append_keyvalue_to_hdf5('resp', respon, elt, output_path + 'resp_', dataset)
        finally:
            lock.release()       

    return domain,respon


def analyse(
        store: Store, input_file: str, output_path: str, task_name: str) -> Store:

    config_tuple = create_configfile(parse_config,input_file, output_path)

    if '' in store.Solution:
        T,M,S,D,R,Neltoverlambda,nRD,g = config_tuple
        dataroot = g.convergemod + '_' + str(g.nff) + '_' + str(int(g.iff*g.model_scale*100)) + '_' + str(int(g.fff*g.model_scale*100)) 
    else:
        config  = store.Configuration
        dataroot = list(config.keys())[0]
        g = config[dataroot]
        T,M,S,D,R,Neltoverlambda,nRD,_ = config_tuple
        

    ppt_per_surface = 1 + int(np.floor(g.piston_radius * (Neltoverlambda / 100)))

    k0, kr, kr_length, dr, dr_length, sfr, RD, lambda0, d_cur = heuristic(nRD, g)


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
    task_q = JoinableQueue()
    resp_q = JoinableQueue()
    # Create processes
    processes = []
    num_proc = min(g.ffu-g.ifu,cpc)
    
    for i in range(num_proc):
        p = Process(target=worker, args=(task_q,resp_q,apod,lock))
        p.start()
        processes.append(p)

    for jj in range(kr_length):
        for pp in range(dr_length):
            path[2] = '_' + str(int(g.skr[jj])) + '_' + str(int(g.sdr[pp]))

            dataset = path[1] + path[2]

            if task_name == "analyse" or task_name == "mfsolution_analyse":
                create_keysized_to_hdf5('doma', domaset_size, output_path + 'doma_', dataset)
                create_keysized_to_hdf5('resp', respset_size, output_path + 'resp_', dataset)
            else:
                store.init_store_behavior('doma_' + dataset,domaset_size)
                store.init_store_behavior('resp_' + dataset,respset_size)

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

                task_q.put((elt,MH,MR,path,task_name,dataset), block=True)

    # Add sentinels to shut down processes.
    print('Adding sentinels...')
    for _ in range(num_proc):
        task_q.put(None)

    for _ in range(kr_length*kr_length*(g.ffu-g.ifu+1)):
        elt,result,dataset = resp_q.get()
        
        if '' not in store.Behavior:
            domain,respon = result
            result = {'doma_'+ dataset: domain, 'resp_'+ dataset: respon}
            store.load_dict_to_store_behaviour(result , elt)
        resp_q.task_done()

    resp_q.join()  # Block until all tasks are done.
    task_q.join() # Block until all elements have been processed.

    return store

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        raise ValueError('Invalid number of arguments. Usage: {} input_file.yaml /path/to/output'.format(sys.argv[0]))

    input_file = sys.argv[1]
    output_path = sys.argv[2]

    store = Store('','')
    analyse(store,input_file, output_path,"analyse")