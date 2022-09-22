import numpy as np
import sys

from scipy.signal import (convolve,convolve2d)
from util.heuristic import heuristic
from util.h5py_util import *
#from threading import Thread
from multiprocessing import Process,Lock,JoinableQueue
from os import cpu_count,getpid
from config import parse_config, create_configfile

def worker(q,l):
    """The process will continually pull elements from the shared queue 
    to process until reaching a None sentinel.
    """

    while True:

        current = q.get(timeout=5)  
        if current is None:
            q.task_done()
            break
        elt,apod,path = current
        print("Starting to process elt: {}".format(elt))
        fn_analyse(elt,apod,l,path)
  
        print("Finished processing elt: {}".format(elt))
        q.task_done()

def fn_analyse(elt,apod,l,path):
    
    output_path,dataroot,heurisset = path

    datafile = dataroot + '_' + str(elt+1) + heurisset
    dataset  = dataroot + heurisset

    apod2 = apod[:,np.newaxis] @ apod[np.newaxis,:]

    with h5py.File(datafile + '.h5', 'r') as f:
        MH = f['domain']
        MR = f['receiver']

        print('DataFile {} read'.format(datafile))

        domain = convolve(MH,apod[np.newaxis,:],mode='valid')
        response = convolve2d(MR,apod2,mode='valid')

    l.acquire()
    try:
        append_keyvalue_to_hdf5('doma', domain, elt, output_path + 'doma_', dataset)
        append_keyvalue_to_hdf5('resp', response, elt, output_path + 'resp_', dataset)
    finally:
        l.release()

    print('DataFile {} appended'.format(output_path + '_' + dataset))


def analyse(config_file, output_path):

    config_tuple = create_configfile(parse_config,config_file)

    T,S,D,R,Neltoverlambda,nRD,g = config_tuple

    ppt_per_surface = 1 + int(np.floor(g.piston_radius * (Neltoverlambda / 100)))

    k0, kr, kr_length, dr, dr_length, sfr, RD, lambda0, d_cur = heuristic(nRD, g)

    dataroot = g.convergemod + '_' + str(g.nff) + '_' + str(int(g.iff*g.model_scale*100)) + '_' + str(int(g.model_scale*100))

    rshape = (R.c.size-ppt_per_surface+1,T.c.size-ppt_per_surface+1)
    dshape = (D.c.size,T.c.size-ppt_per_surface+1)

    domaset_size = (sfr.size,) + dshape
    respset_size = (sfr.size,) + rshape

    range_ppt = .5+np.arange(-int(ppt_per_surface/2),int(ppt_per_surface/2))

    apod = 8/(ppt_per_surface*np.pi)*np.sqrt(int(ppt_per_surface/2)**2-range_ppt**2)
    
    processes = []

    cpc = cpu_count()

    lock = Lock()
    
    path = [output_path,dataroot,'_']

    queue_max_size = 20
    q = JoinableQueue(maxsize=queue_max_size)
    # Create processes
    processes = []
    
    for i in range(cpc):
        p = Process(target=worker, args=(q,lock))
        p.start()
        processes.append(p)

    for jj in range(kr_length):
               
        for pp in range(dr_length):
            
            path[2] = '_' + str(int(g.skr[jj])) + '_' + str(int(g.sdr[pp]))

            dataset = path[1] + path[2]

            create_keysized_to_hdf5('doma', domaset_size, output_path + 'doma_', dataset)
            create_keysized_to_hdf5('resp', respset_size, output_path + 'resp_', dataset)

            for elt in range(g.ifu-1,g.ffu):

                #fn_analyse(config_tuple,datafile,elt,dataset)
                    # create processes
                #processes.append(Process(target=fn_analyse, args=(elt,apod,lock,path)))
                #q.put((elt,apod,lock,path), block=True)
                q.put((elt,apod,path), block=True)

    # Add sentinels to shut down processes.
    print('Adding sentinels...')
    for _ in range(cpc):
        q.put(None)

    q.join() # Block until all elements have been processed.

if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        raise ValueError('Invalid number of arguments. Usage: {} input_file.yaml /path/to/output'.format(sys.argv[0]))

    input_file = sys.argv[1]
    output_path = sys.argv[2]

    analyse(input_file, output_path)