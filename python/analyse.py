import numpy as np
import sys

from scipy.signal import (convolve,convolve2d)
from util.heuristic import heuristic
from util.h5py_util import *
from config import *
import time


def fn_analyse(config_tuple,datafile):
    
    T,_,D,R,Neltoverlambda,nRD,g = config_tuple
    
    response = structtype()
    field = structtype()
    
    ppt_per_surface = 1 + int(np.floor(g.piston_radius * (Neltoverlambda / 100)))

    range_ppt = .5+np.arange(-int(ppt_per_surface/2),int(ppt_per_surface/2))
    apodization = 8/(ppt_per_surface*np.pi)*np.sqrt(int(ppt_per_surface/2)**2-range_ppt**2)

    apodization2 = apodization[:,np.newaxis] @ apodization[np.newaxis,:]
    '''
    range_ppt = np.arange(ppt_per_surface)
    
    ndx = structtype()
    ndx.T = structtype()
    ndx.R = structtype()
      
    response.range = structtype()
    response.pool = structtype()
    response.pitch = structtype()
    response.catch = structtype()
    response.ndx = structtype()
    
    response.pitch.size = (T.a.size - ppt_per_surface + 1)
    response.catch.size = (R.a.size - ppt_per_surface + 1) 
    response.pitch.A = np.zeros(response.pitch.size)
    response.catch.A = np.zeros(response.catch.size)

    #response.pool.p = np.zeros(response.catch.size, dtype=np.complex128)    
    response.p = np.zeros((response.pitch.size,response.catch.size), dtype=np.complex128)    
    field.p = np.zeros((response.pitch.size,D.c.size), dtype=np.complex128)    
    '''

    with h5py.File(datafile + '.h5', 'r') as f:
        MH = f['domain']
        MR = f['receiver']

        response.pitch = convolve(T.z,apodization,mode='valid')
        response.catch = convolve(R.z,apodization,mode='valid')

        field.p = convolve(MH,apodization[np.newaxis,:],mode='valid')
        response.p = convolve2d(MR,apodization2,mode='valid')
        '''
        for t in np.arange(0,response.pitch.size):
            
            #   EXTRACT TRANSMITTER PISTON INDEXES          
            t_range = t + range_ppt
            response.pitch.A[t] = np.mean(Tz[t_range])

            #   FILTER FIELD PARAMETERS   
            field.p[t,...] = np.sum(MH[...,t_range], 1)
            response.pool.p = np.sum(MR[...,t_range], 1)
                        
            for r in np.arange(0,response.catch.size):

                #   EXTRACT RECEIVER PISTON INDEXES
                r_range = r + range_ppt

                response.catch.A[r] = np.mean(Rz[r_range])
   
                #   ANALYSE RESPONSE RANGES                  
                response.p[t,r] = np.sum(response.pool.p[r_range])

        '''
        resp = np.copy(response.p)
        doma = np.copy(field.p)
        
    return  resp,doma


def analyse(fn_analyse, config_file, output_path):
    config_tuple = create_configfile(parse_config,config_file)

    T,S,D,R,_,nRD,g = config_tuple

    k0, kr, kr_length, dr, dr_length, sfr, RD, lambda0, d_cur = heuristic(nRD, g)

    dataroot = g.convergemod + '_' + str(g.nff) + '_' + str(int(g.iff*g.model_scale*100)) + '_' + str(int(g.model_scale*100))

    for jj in range(kr_length):
               
        for pp in range(dr_length):
            
            dataset = dataroot + '_' + str(int(g.skr[jj])) + '_' + str(int(g.sdr[pp]))
            
            for ii in range(g.ifu-1,g.ffu):
                
                datafile = dataroot + '_' + str(ii+1) + '_' + str(int(g.skr[jj])) + '_' + str(int(g.sdr[pp]))

                stop_time = time.time()
                start_time = stop_time

                r,d = fn_analyse(config_tuple,datafile)

                stop_time = time.time()
                print(stop_time - start_time)
                if ii==0:
                    domaset_size = (sfr.size,) + d.shape
                    respset_size = (sfr.size,) + r.shape
                    create_keysized_to_hdf5('doma', domaset_size, output_path + 'doma_', dataset)
                    create_keysized_to_hdf5('resp', respset_size, output_path + 'resp_', dataset)           
                append_keyvalue_to_hdf5('doma', d, ii, output_path + 'doma_', dataset)
                append_keyvalue_to_hdf5('resp', r, ii, output_path + 'resp_', dataset)


if __name__ == "__main__":
    
    if len(sys.argv) != 3:
        raise ValueError('Invalid number of arguments. Usage: {} input_file.yaml /path/to/output'.format(sys.argv[0]))

    input_file = sys.argv[1]
    output_path = sys.argv[2]

    analyse(fn_analyse, input_file, output_path)