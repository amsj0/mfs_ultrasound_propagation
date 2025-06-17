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


def mfsolution(name,config_file, output_path):

    print(config_file)
    print(output_path)
    T,P,S,nRD,g = reconfigure(create_configfile(parse_config,config_file, output_path))  
         
    CP = Compute(P,T,S)

    CP.InitCL("GPU")

    k0, kr, kr_length, dr, dr_length, sfr, RD, lambda0, d_cur = heuristic(nRD, g)
   
    dataroot = g.convergemod + '_' + str(g.nff) + '_' + str(int(g.iff*g.model_scale*100)) + '_' + str(int(g.fff*g.model_scale*100)) 

    datafile = ''
    
    ST = Store(datafile)

    l = 2

    for ii in range(g.ifu-1,g.ffu):
        
        k_cur = k0[ii]*lambda0/(RD)
        
        sf_cur = sfr[ii]
        print('Spectrum Ratio')
        print(sf_cur)


        k_out = k_cur
        d_out = d_cur
        p_out = [k_out,d_out]
        CP.compute_lower_side(p_out)

        for jj in range(kr_length):
            
            k_r = kr[jj]

            for pp in range(dr_length):
                
                d_aux = d_out*dr[pp]

                d_r = dr[pp]
                dc_cur = g.rjR(k_r,d_r)[ii]
                kc_cur = g.keq(k_r,d_r)[ii]
                d_r = dc_cur/d_out

                            
                p_cur = [kc_cur,dc_cur]

                datafile = dataroot + '_' + str(ii+1) + '_' + str(int(g.skr[jj])) + '_' + str(int(g.sdr[pp]))

                print('Wavenumber Ratio')
                print(kr[jj]/(1-1j*g.att)/np.abs(1-0*1j*g.att))
                print('Density Ratio')
                print(d_r)
                
                CP.compute_upper_side(p_cur,p_out)
                CP.propagate_transfer()

                if name == "__main__":
                    ST.load_dict_to_store(CP.M, datafile)
                elif name == "mfsolution":
                    save_dict_to_hdf5(CP.M, output_path, datafile)
                
                print('Datafile {} created'.format(datafile))
    return ST


if __name__ == "__main__":

    if len(sys.argv) != 3:
        raise ValueError('Invalid number of arguments. Usage: {} input_file.yaml'.format(sys.argv[0]))
    
    input_file = sys.argv[1]
    output_path = sys.argv[2]    

    mfsolution("mfsolution",input_file, output_path)