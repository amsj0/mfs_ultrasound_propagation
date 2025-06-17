#from nbformat import read
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
from scipy.signal import hilbert
import scipy.signal as sp
from util.spectrum import Spectrum
import yaml

import sys

from util.h5py_util import *

def load_dataset(conve_mod, gridname, pathname, filename):
    
    grid,respc,scale,ndx0 = load_para(pathname, conve_mod, gridname)
    
    data_set = {
        'doma' : load_(pathname,'doma',conve_mod + filename),
        'resp' : load_(pathname,'resp',conve_mod + filename)
    }

    return data_set,grid,respc,scale,ndx0


def pre_config(config_file,output_path):
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

        spec_band = cfg['spec_band']
        spec_size = cfg['spec_size']
        ref_freq = cfg['ref_freq']
        conve_mod = cfg['CONVE_MOD']
        
    gridname = '_' + str(numbr_freq) + '_' + str(initi_freq) + '_' + str(final_freq) 
    
    filename = gridname + '_' + str(int(skr)) + '_' + str(int(sdr)) + '.h5'    
              
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

    expr = np.loadtxt(output_path + 'experimental_pulse.txt', dtype=np.complex_)
    
    f_dom = np.linspace(1/x_size,1,x_size) - 0.5

    central_range = ref_freq*(1*f_dom + final_freq/(100*2))

    data_set,grid,respc,scale,ndx0 = load_dataset(conve_mod, gridname, output_path, filename)

    SP = Spectrum(initi_freq,final_freq,numbr_freq,parti_freq,x,filter,expr,spec_size)

    return SP, offset, dtsr, x_size, central_range, data_set, grid,respc,scale , ndx0, x_spec_full, conve_mod

def plt_arrange(i):
    
    _,_,tseries = SP.synth_fseries_from_centr_freq(central_range[i])

    plt.plot(np.max(abs(tseries),axis=0).transpose())
    plt.savefig('tseries' + str(i) + '.png')
    plt.clear(True)


def plt_max_tseries(pathname,converge):
    
    resp = load_(pathname,'resp',converge)

    new_resp = SP.expand_resp_2(resp)
    max_tseries = np.empty((x_size,new_resp.shape[2],new_resp.shape[1]),dtype='double')
     
    for i in range(x_size):
        _,_,tseries = SP.synth_tseries(new_resp,central_range[i])
        max_tseries[i,:,:] = np.max(abs(tseries),axis=0).transpose()

    return max_tseries


def plt_mat_tseries_1(pathname,converge):
    
    resp = load_(pathname,'resp',converge)
    
    mat_tseries = np.empty((x_size,resp.shape[2]),dtype='complex')
     
    for i in range(int(x_size)):
        _,_,spec = SP.synth_fseries_from_centr_freq(central_range[i])
        for j  in range(int(resp.shape[2])):
            new_resp = SP.expand_resp_1(resp[:,:,j])
            _,tseries = SP.synth_tseries_from_spec(new_resp,spec)
            tsumseries = np.sum(tseries,axis=-1)
            tstdseries = np.std(tsumseries)
            mat_tseries[i,j] = tstdseries
    return mat_tseries

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

def create_matrix(SP, dtsr, x_size, central_range, grid, data_set, ndx0, x_spec_full, mat_tseries, offset, vlim, rlim, vlimmax):
    respttm = []

    if dtsr:
        this_series, fg, h0 = set_domain_plot(grid,ndx0,data_set)


    for i in range(int(1*x_size*4/8-0)-1,int(1*x_size*4/8-0),1):
        central_freq = central_range[i]
        osc,freq,spec = SP.synth_fseries_from_centr_freq(central_freq)
        #osc,freq,spec = SP.synth_fseries_from_experiment()
        spec0 = spec[0:int(SP.spec_size)]
        
        for j in range(data_set['resp'].shape[-1]):
            resp_data = data_set['resp'][...,j]
            
            for k in range(data_set['resp'].shape[-2]):
                resptt = resp_data[...,k].transpose()
                
                if k==j:
                      if (k==0):
                        respttm.append(resptt)


                new_resp = np.conj(SP.expand_resp_new(resptt/respttm[0]))
                # new_resp = np.conj(SP.expand_resp_new(resptt))
                tseries = SP.synth_tseries_from_spec_full_new(new_resp*spec0)    
                # tseries = SP.synth_tseries_from_spec_full_filtered(new_resp*spec0)    
        
                vec_tseries = tseries.transpose()


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
    

def plot_data(x_spec_full, offset, vlim, rlim, vlimmax, probv):
    
    spec_grid,resp_grid = np.meshgrid(x_spec_full,probv)
    tx_grid,rx_grid = np.meshgrid(probv,probv)

    x_grid = (tx_grid*np.cos(np.pi/4) - rx_grid*np.sin(np.pi/4))/np.sqrt(2)
    y_grid = (tx_grid*np.sin(np.pi/4) + rx_grid*np.cos(np.pi/4))/np.sqrt(2)

    fg,axs = plt.subplots(2,1,figsize=(10,7))
 
    axs[0].set_title('Diagonal Signal')
    axs[0].pcolormesh(spec_grid,resp_grid,np.abs(vlim),shading='nearest')

    vlimmaxN = vlimmax/np.max(vlimmax[0])

    o_grid = []
    v_grid = []

    for off in np.arange(-int(probv.size/2),1+int(probv.size/2),5):
        o_grid.append(np.diagonal(x_grid,offset=off))
        v_grid.append(np.diagonal(y_grid,offset=off))
        v_grid.append(np.diagonal(vlimmaxN,offset=off))


    t_grid = np.diagonal(y_grid,offset=offset)
    n_grid = np.diagonal(x_grid,offset=offset)
    tpec_grid,tesp_grid = np.meshgrid(x_spec_full,t_grid)
    
    axs[1].set_title('{}-Off Diagonal Signal'.format(probv[int(probv.size/2)+offset]))
    axs[1].pcolormesh(tpec_grid,tesp_grid,np.abs(rlim),shading='nearest')

    plt.figure(figsize=(7,7))
    plt.pcolormesh(tx_grid,rx_grid,vlimmaxN,shading='nearest',cmap='Greys')

    x_table = probv
    y_table = np.diagonal(vlimmaxN)
    y_off_table = np.diagonal(vlimmaxN,offset=offset)


    plt.figure('diagonal_matrix')
    plt.plot(x_table,y_table)
    plt.plot(t_grid,y_off_table)
    plt.ylim(0,1)
    plt.title('Off Diagonal Matrix')
    plt.tight_layout()

    plt.figure('full_off_diagonal_matrix')
    plt.plot(*v_grid)
    plt.title('Full Off Diagonal Matrix')
    plt.tight_layout()    

def save_table(conve_mod, x_spec_full, freq, vlim, vlimmax, probv, offset):
    
    vlimmaxN = vlimmax/np.max(vlimmax[0]) # Normalize the maximum value of the diagonal signal
    #vlimmaxN = vlimmax/np.max(vlimmax[-1])

    tx_grid,rx_grid = np.meshgrid(probv,probv)

    y_grid = (tx_grid*np.sin(np.pi/4) + rx_grid*np.cos(np.pi/4))/np.sqrt(2)

    t_grid = np.diagonal(y_grid,offset=offset)
    y_off_table = np.diagonal(vlimmaxN,offset=offset)

    d_table = np.array([probv,np.diagonal(vlimmaxN)]).transpose()
    d_off_table = np.array([t_grid,y_off_table]).transpose()
    d_time_table = np.block([[0,x_spec_full.transpose()],[probv[:,np.newaxis],np.real(vlim)]])

    np.savetxt('T'+conve_mod+'_'+str(freq)+'.csv',d_time_table, delimiter=',')
    np.savetxt('V'+conve_mod+'_'+str(freq)+'.csv',d_table, delimiter=',', header=','.join(('height','amax')), comments='')
    np.savetxt('M'+conve_mod+'_'+str(freq)+'.csv',vlimmaxN, delimiter=',')
    np.savetxt('Voff'+conve_mod+'_'+str(freq)+'.csv',d_off_table, delimiter=',', header=','.join(('height','amax')), comments='')


def set_empty_matrix(SP, offset, data_set, respc, scale):
    mat_tseries = np.empty((SP.spec_size,data_set['doma'].shape[-2]),dtype='complex')
        
    vlim = np.empty((data_set['resp'].shape[-1],SP.spec_size*2),dtype='complex')
    rlim = np.empty(((data_set['resp'].shape[-1]-np.abs(offset)),SP.spec_size*2),dtype='complex')
    
    vlimmax = np.empty((int(data_set['resp'].shape[-2]),int(data_set['resp'].shape[-1])),dtype='float')

    probv =  respc+(.5+np.arange(-int(data_set['resp'].shape[-1])/2,int(data_set['resp'].shape[-1])/2))[:]*scale
    return mat_tseries,vlim,rlim,vlimmax,probv

def tseries(name,config_file, output_path):

    SP, offset, dtsr, x_size, central_range, data_set, grid ,respc ,scale , ndx0, x_spec_full, conve_mod = pre_config(config_file,output_path)

    mat_tseries, vlim,rlim, vlimmax, probv = set_empty_matrix(SP, offset, data_set, respc, scale)

    freq = create_matrix(SP, dtsr, x_size, central_range, grid, data_set, ndx0, x_spec_full, mat_tseries, offset, vlim,rlim, vlimmax)

    save_table(conve_mod, x_spec_full, freq, vlim, vlimmax, probv, offset)

    #plot_data(x_spec_full, vlim, vlimmax, probv)

    plt.show()


if __name__ == '__main__':
   
    if len(sys.argv) != 3:
        raise ValueError('Invalid number of arguments. Usage: {} /path/to/output  config.yaml'.format(sys.argv[0]))

    output_path = sys.argv[1]
    config_file = sys.argv[2]

    tseries("tseries",config_file, output_path)