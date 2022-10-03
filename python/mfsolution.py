import numpy as np
import sys

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


def mfsolution(name,config_file):

    T,P,S,nRD,g = reconfigure(create_configfile(parse_config,config_file))  
         
    CP = Compute(P,T,S)

    CP.InitCL("GPU")

    k0, kr, kr_length, dr, dr_length, sfr, RD, lambda0, d_cur = heuristic(nRD, g)
   
    dataroot = g.convergemod + '_' + str(g.nff) + '_' + str(int(g.iff*g.model_scale*100)) + '_' + str(int(g.model_scale*100)) 

    datafile = ''
    
    ST = Store(datafile)

    l = 2

    for ii in range(g.ifu-1,g.ffu):
        
        k_cur = k0[ii]*lambda0/(RD)
        
        sf_cur = sfr[ii]
        print('Spectrum Ratio')
        print(sf_cur)

        CP.compute_lower_side(k_cur)

        k_out = k_cur
        d_out = d_cur

        for jj in range(kr_length):
            
            k_r = kr[jj]
            k_curi = k_out/k_r
            k_cur  = k_curi+1j*g.att*g.sfr[ii]**(1+l)*g.sfr[int(g.ffu/2)]**(1-l)

            for pp in range(dr_length):
                
                d_r = dr[pp]
                d_cur  = d_out/dr[pp]
                            
                datafile = dataroot + '_' + str(ii+1) + '_' + str(int(g.skr[jj])) + '_' + str(int(g.sdr[pp]))

                print('Wavenumber Ratio')
                print(kr[jj]/(1-1j*g.att)/np.abs(1-0*1j*g.att))
                print('Density Ratio')
                print(dr[pp])
                
                CP.compute_upper_side(k_cur,k_out,d_r)
                CP.propagate_transfer()

                if name == "__main__":
                    ST.load_dict_to_store(CP.M, datafile)
                elif name == "mfsolution":
                    save_dict_to_hdf5(CP.M, datafile)
                #save_keyvalue_to_hdf5('doma',CP.M['domain'], '', datafile)
                
                print('Datafile {} created'.format(datafile))
    return ST


if __name__ == "__main__":

    if len(sys.argv) != 2:
        raise ValueError('Invalid number of arguments. Usage: {} input_file.yaml'.format(sys.argv[0]))
    
    input_file = sys.argv[1]

    mfsolution("mfsolution",input_file)