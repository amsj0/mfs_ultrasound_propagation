import numpy as np
import sys

from scipy.signal import (convolve,convolve2d)
from util.heuristic import heuristic
from util.h5py_util import *
from config import *
import time


def fn_analyse(config_tuple,datafile):
    
    T,_,D,R,Neltoverlambda,nRD,g = config_tuple
    #loadmat([path,datafile])
    #loadmat([path,configfile + '.h5'])
    #loadmat([path,analysisfile])

    #   DEFINES PISTON INDEXES
    #Tz = T.z
    #Rz = R.z
    
    response = structtype()
    field = structtype()
    
    ppt_per_surface = 1 + int(np.floor(g.piston_radius * (Neltoverlambda / 100)))

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

        response.pitch = convolve(T.z,np.ones(ppt_per_surface),mode='same')
        response.catch = convolve(R.z,np.ones(ppt_per_surface),mode='same')

        field.p = convolve(MH,np.ones((1,ppt_per_surface)),mode='same')
        response.p = convolve2d(MR,np.ones((ppt_per_surface,ppt_per_surface)),mode='same')
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


def analyse(fn_analyse, yaml_path):
    config_tuple = create_configfile(parse_config,yaml_path)

    T,S,D,R,_,nRD,g = config_tuple

    k0, kr, kr_length, dr, dr_length, sfr, RD, lambda0, d_cur = heuristic(nRD, g)

    resp = []
    doma = []
    path = "D:/MATLAB/menisco/"

    analysisfile = "P.mat"

    dataroot = g.convergemod + '_' + str(g.nfr) + '_' + str(int(g.iff*g.model_scale*100)) + '_' + str(int(g.model_scale*100))

    for jj in range(kr_length):
               
        for pp in range(dr_length):
            
            for ii in range(g.nfi-1,g.nff):
                
                datafile = dataroot + '_' + str(ii+1) + '_' + str(int(g.skr[jj])) + '_' + str(int(g.sdr[pp]))

                stop_time = time.time()
                start_time = stop_time

                r,d = fn_analyse(config_tuple,datafile)
                resp.append(r)
                doma.append(d)

                stop_time = time.time()
                print(stop_time - start_time)
                
            dataset = dataroot + '_' + str(int(g.skr[jj])) + '_' + str(int(g.sdr[pp]))

            save_keyvalue_to_hdf5('doma', doma, path + 'doma_', dataset)
            save_keyvalue_to_hdf5('resp', resp, path + 'reps_', dataset)

            resp = []
            doma = []


if __name__ == "__main__":
    
    if len(sys.argv) != 2:
        raise ValueError('Invalid number of arguments. Usage: {} config_file.yaml'.format(sys.argv[0]))

    yaml_path = sys.argv[1]

    analyse(fn_analyse, yaml_path)