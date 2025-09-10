#from nbformat import read
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
from scipy.signal import hilbert
import scipy.signal as sp
from multiprocessing import Process, JoinableQueue, cpu_count
from util.spectrum import Spectrum
from util.store import Store
from util.entry import Entry
from util.h5py_util import *
import yaml

import sys

from util.h5py_util import *

def load_dataset(task_name,store ,conve_mod, gridname, pathname, filename):
    
    grid,respc,scale,ndx0 = load_para(pathname, conve_mod, gridname)
    
    dataset = conve_mod + filename
    if task_name == "tseries":
        data_set = {
            'doma' : load_(pathname,'doma',dataset + '.h5'),
            'resp' : load_(pathname,'resp',dataset + '.h5')
        }
    else:
        data_set = {
            'doma' : store.Behavior['doma_' + dataset],
            'resp' : store.Behavior['resp_' + dataset]
        }

    return data_set,grid,respc,scale,ndx0


def pre_config(task_name, store, config_file,output_path):
    """
    Prepares the configuration for the ultrasound propagation simulation.
    Parameters:
        config_file (str): Path to the configuration file in YAML format.
        output_path (str): Directory path where output files will be saved.
    Returns:
        tuple: A tuple containing:
            - SP (Spectrum): An instance of the Spectrum class initialized with simulation parameters.
            - offset (float): Offset value from the configuration.
            - dtsr (float): DTSR value from the configuration.
            - x_size (int): Size of the x array.
            - central_range (numpy.ndarray): Central frequency range calculated from reference frequency.
            - data_set (numpy.ndarray): Loaded dataset from the specified model.
            - grid (numpy.ndarray): Grid data from the loaded dataset.
            - respc (numpy.ndarray): Response data from the loaded dataset.
            - scale (float): Scale factor from the loaded dataset.
            - ndx0 (int): Index value from the loaded dataset.
            - x_spec_full (numpy.ndarray): Full spectrum x array.
            - conve_mod (str): Conversion model from the configuration.
            - tranx_num (int): Number of transducers from the configuration.
    Note:
        freq_scale is calculated as the ratio of the specified spectrum size to the number of frequencies divided by the spectrum band.
    """
   
    with open(config_file, 'r') as f:
        cfg = yaml.safe_load(f)
        initi_freq = cfg['initi_freq']
        final_freq = cfg['final_freq']
        parti_freq = cfg['parti_freq']
        numbr_freq = cfg['numbr_freq']

        skr = cfg['skr']
        sdr = cfg['sdr']
        dtsr = cfg['dtsr']
        offset = cfg['offset']

        tranx_num = cfg['tranx_num']
        spec_band = cfg['spec_band']
        spec_size = cfg['spec_size']
        ref_freq = cfg['ref_freq']
        conve_mod = cfg['CONVE_MOD']
        input_mod = cfg['INPUT_MOD']
        
    gridname = '_' + str(numbr_freq) + '_' + str(initi_freq) + '_' + str(final_freq) 
    
    filename = gridname + '_' + str(int(skr)) + '_' + str(int(sdr))    
              
    mu, sigma = 0, 0.055 # .125 # 0.055 # -.25 .07 (1 MHz), 0 .13 (1.5 MHz)
    
    freq_scale = int(spec_size/(numbr_freq/spec_band))
    sampling_freq = freq_scale*10*ref_freq
    
    x = np.linspace(1/sampling_freq,freq_scale*numbr_freq/sampling_freq,freq_scale*numbr_freq)[...,np.newaxis]
    x_size = np.size(x,axis=0)
    
    x_spec_full = np.linspace(1/sampling_freq,spec_size/sampling_freq,spec_size)[...,np.newaxis]

    echo_size = 4*ref_freq/sampling_freq

    x_dom = np.linspace(1/x_size,1/freq_scale,x_size)

    gauss_filter = 1/(sigma * np.sqrt(2.0 * np.pi)) * np.exp( - (x_dom - mu)**2 / (2 * sigma**2) )[...,np.newaxis]    
    
    cosine_filter =  ((0.5 - 0.5 * np.cos(2.0 * np.pi * x_dom / echo_size)) * (x_dom<echo_size)*(x_dom>0))[...,np.newaxis]

    filter = cosine_filter

    if input_mod == 'experiment':
        expr = np.loadtxt(output_path + 'experimental_pulse.txt', dtype=np.complex_)
    elif input_mod == 'synthetic':
        expr = np.sin(2 * np.pi * ref_freq * x_dom) * gauss_filter  
    
    f_dom = np.linspace(1/x_size,1,x_size) - 0.5

    central_range = ref_freq*(1*f_dom + final_freq/(100*2))

    data_set,grid,respc,scale,ndx0 = load_dataset(task_name, store, conve_mod, gridname, output_path, filename)

    SP = Spectrum(initi_freq,final_freq,numbr_freq,parti_freq,x,filter,expr,spec_size)

    return SP, offset, dtsr, x_size, central_range, data_set, grid,respc,scale , ndx0, x_spec_full, conve_mod, tranx_num

def set_domain_plot(grid,ndx0, data_set):
    
    this_series = np.empty(np.prod(grid.shape))
    this_series[ndx0] = np.nan
    this_series[np.logical_not(ndx0)] = np.zeros(data_set['doma'].shape[1])
    
    this_grid = grid

    fg,axs = plt.subplots(1,1,figsize=(10,7))
    h0 = axs.pcolormesh(this_grid.real,this_grid.imag,this_series.reshape(this_grid.shape),shading='nearest',vmin=-3,vmax=3)
    axs.set_aspect('equal', 'box')
    axs.tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=True)
    axs.set_title('Total Field in Domain')

    fg.canvas.draw()
    fg.canvas.flush_events()

    return this_series, fg, h0

def create_matrix(SP, dtsr, x_size, central_range, grid, data_set, ndx0, x_spec_full, mat_tseries, offset, lim, vlimmax):
    respttm = []

    vlim,mlim,rlim = lim

    if dtsr:
        this_series, fg, h0 = set_domain_plot(grid,ndx0,data_set)

    data_x_size = data_set['resp'].shape[-1]
    data_y_size = data_set['resp'].shape[-2]

    matrix_x_strafe = (data_x_size-1)/(mlim.shape[0])
    matrix_y_strafe = (data_y_size-1)/(mlim.shape[1])

    for i in range(int(1*x_size*4/8-0)-1,int(1*x_size*4/8-0),1):
        central_freq = central_range[i]
        osc,freq,spec = SP.synth_fseries_from_centr_freq(central_freq)
        #osc,freq,spec = SP.synth_fseries_from_experiment()
        spec0 = spec[0:int(SP.spec_size)]
        

        print('Central frequency: ' + str(central_freq) + ' MHz' )
        for j in range(data_x_size):
            resp_data = data_set['resp'][...,j]

            matrix_x_ndx = j % matrix_x_strafe
            
            for k in range(data_y_size):
                resptt = resp_data[...,k].transpose()
                
                matrix_y_ndx = k % matrix_y_strafe

                if k==j:
                      if (k==0):
                        respttm.append(resptt)


                # new_resp = np.conj(SP.expand_resp_new(resptt/respttm[0]))
                new_resp = np.conj(SP.expand_resp_new(resptt))
                tseries = SP.synth_tseries_from_spec_full_new(new_resp*spec0)    
                # tseries = SP.synth_tseries_from_spec_full_filtered(new_resp*spec0)    
        
                vec_tseries = tseries.transpose()

                if not matrix_x_ndx and not matrix_y_ndx:
                    mlim[int(-.5+j/matrix_x_strafe),int(-.5+k/matrix_y_strafe),...] = vec_tseries

                if k==j:
                    vlim[j,...] = vec_tseries
                if k==(j+offset):
                    rlim[j - int((np.abs(offset) - offset)/2),...] = vec_tseries                    

                vlimmax[k,j] = np.max(np.abs(vec_tseries),axis=-1)
                


            if dtsr:
                doma_data = data_set['doma'][...,j]
                            
                for k in range(int(doma_data.shape[-1])):
                    domatt = doma_data[...,k].transpose()
                    new_doma = SP.expand_resp_new(domatt)
                    tseries = SP.synth_tseries_from_spec_full_new(new_doma*spec0)
                    mat_tseries[...,k] = tseries.transpose()

                mlim = np.amax(np.abs(mat_tseries))/4

                for ll in range(int(SP.spec_size*2/dtsr)):
                    this_series[np.logical_not(ndx0)] = mat_tseries[(ll+1)*dtsr-1,:].real/mlim

                    h0.set_array(this_series.ravel())
                    
                    fg.suptitle('time-step ' + str((ll+1)*dtsr-1) + '| time ' + str(x_spec_full[(ll+1)*dtsr-1]) + ' | height-step ' + str(j))
                    fg.canvas.draw()
                    fg.canvas.flush_events()
    return int(central_freq/central_range[int(x_size*1/2)-1]*100)
    

def save_table(conve_mod, x_spec_full, freq, lim, vlimmax, prob, offset, output_path=''):
    
    vlimmaxN = vlimmax/np.max(vlimmax[0]) # Normalize the maximum value of the diagonal signal
    #vlimmaxN = vlimmax/np.max(vlimmax[-1])

    vlim,mlim,rlim = lim
    probv,probm = prob

    tx_grid,rx_grid = np.meshgrid(probv,probv)

    y_grid = (tx_grid*np.sin(np.pi/4) + rx_grid*np.cos(np.pi/4))/np.sqrt(2)

    t_grid = np.diagonal(y_grid,offset=offset)
    y_off_table = np.diagonal(vlimmaxN,offset=offset)

    d_table = np.array([probv,np.diagonal(vlimmaxN)]).transpose()
    d_off_table = np.array([t_grid,y_off_table]).transpose()
    d_time_table = np.block([[0,x_spec_full.transpose()],[probv[:,np.newaxis],np.real(vlim)]])
    
    m_time_table = np.block([[0,0,x_spec_full.transpose()],[probm.reshape(probm.shape[0]*probm.shape[1],probm.shape[2]),np.real(mlim.reshape(mlim.shape[0]*mlim.shape[1],mlim.shape[2]))]])


    data = {
        "T": Entry(d_time_table, ("height", "a")),
        "V": Entry(d_table, ("height", "amax")),
        "N": Entry(m_time_table, ("height","height", "a")),
        "M": Entry(vlimmaxN, ("height", "amax")),
        "Voff": Entry(d_off_table, ("height", "amax"))
    }
    
    save_data_to_hdf5(
        data,
        output_path,
        {'conve_mod': conve_mod,'freq': freq}
        )
    np.savetxt(output_path+'T'+conve_mod+'_'+str(freq)+'.csv',d_time_table, delimiter=',')
    np.savetxt(output_path+'N'+conve_mod+'_'+str(freq)+'.csv',m_time_table, delimiter=',')
    np.savetxt(output_path+'V'+conve_mod+'_'+str(freq)+'.csv',d_table, delimiter=',', header=','.join(('height','amax')), comments='')
    np.savetxt(output_path+'M'+conve_mod+'_'+str(freq)+'.csv',vlimmaxN, delimiter=',')
    np.savetxt(output_path+'Voff'+conve_mod+'_'+str(freq)+'.csv',d_off_table, delimiter=',', header=','.join(('height','amax')), comments='')


def set_empty_matrix(SP, offset, data_set, respc, scale, tranx_num=1):
    mat_tseries = np.empty((SP.spec_size,data_set['doma'].shape[-2]),dtype='complex')
        
    data_m_size = data_set['resp'].shape[-2],data_set['resp'].shape[-1]
    respv_size = (np.array(data_m_size) - 1) / 2

    probmn = tuple()
    probvn = tuple()
    for respv_s in respv_size:
        respvn = 0 * respc + (np.arange(-respv_s, 1 + respv_s))[:] * scale
        step_n = int(2 * respv_s / tranx_num)
        start_n = len(respvn) // 2 - (tranx_num // 2) * step_n + int(step_n / 2)
        stop_n = start_n + tranx_num * step_n
        probmn += (respvn[start_n:stop_n:step_n],)
        probvn += (respvn,)

    
    # respv_x = 0 * respc + (np.arange(-respv_size[1], 1 + respv_size[1]))[:] * scale
    # step_x = int(2 * respv_size[1] / tranx_num)
    # start_x = len(respv_x) // 2 - (tranx_num // 2) * step_x
    # stop_x = start_x + tranx_num * step_x
    # probmx = respv_x[start_x:stop_x:step_x]

    # respv_y = 0 * respc + (np.arange(-respv_size[0], 1 + respv_size[0]))[:] * scale
    # step_y = int(2 * respv_size[0] / tranx_num)
    # start_y = len(respv_y) // 2 - (tranx_num // 2) * step_y
    # stop_y = start_y + tranx_num * step_y
    # probmy = respv_y[start_y:stop_y:step_y]
    probt = np.meshgrid(*probmn)
    probm = np.stack(probt,axis=-1)
    probv = probvn[0]

    vlim = np.empty((probv.shape[0],SP.spec_size*2),dtype='complex')
    mlim = np.empty((probm.shape[0],probm.shape[1],SP.spec_size*2),dtype='complex')
    rlim = np.empty(((data_set['resp'].shape[-1]-np.abs(offset)),SP.spec_size*2),dtype='complex')
    
    vlimmax = np.empty((int(data_set['resp'].shape[-2]),int(data_set['resp'].shape[-1])),dtype='float')
    
    return mat_tseries,(vlim,mlim,rlim),vlimmax,(probv,probm)

def tseries(store,config_file, output_path,task_name):

    SP, offset, dtsr, x_size, central_range, data_set, grid ,respc ,scale , ndx0, x_spec_full, conve_mod, tranx_num = pre_config(task_name, store, config_file,output_path)

    mat_tseries, lim, vlimmax, prob = set_empty_matrix(SP, offset, data_set, respc, scale, tranx_num)

    freq = create_matrix(SP, dtsr, x_size, central_range, grid, data_set, ndx0, x_spec_full, mat_tseries, offset, lim, vlimmax)

    save_table(conve_mod, x_spec_full, freq, lim, vlimmax, prob, offset, output_path)

if __name__ == '__main__':
   
    if len(sys.argv) != 3:
        raise ValueError('Invalid number of arguments. Usage: {} /path/to/output  config.yaml'.format(sys.argv[0]))

    output_path = sys.argv[1]
    config_file = sys.argv[2]


    store = Store('','')
    tseries(store,config_file, output_path,"tseries")