#from nbformat import read
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
#from scipy.signal import hilbert
from util.spectrum import Spectrum
import yaml

import sys

from util.h5py_util import *

def load_dataset(conve_mod, gridname, pathname, filename):
    
    grid,respc,scale,ndx0 = load_para('', conve_mod, gridname)
    
    data_set = {
        'doma' : load_(pathname,'doma',conve_mod + filename),
        'resp' : load_(pathname,'resp',conve_mod + filename)
    }

    return data_set,grid,respc,scale,ndx0


def pre_config(config_file,output_path):
   
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

        spec_size = cfg['spec_size']
        ref_freq = cfg['ref_freq']
        conve_mod = cfg['CONVE_MOD']
        radia_siz = cfg['RADIA_SIZ']
        eleme_wav = cfg['ELEME_WAV']
       
    gridname = '_' + str(numbr_freq) + '_' + str(initi_freq) + '_' + str(final_freq) 
    
    filename = gridname + '_' + str(skr) + '_' + str(sdr) + '.h5'    
           
    grid_scale = radia_siz*eleme_wav/1000
    
    mu, sigma = 0, 0.05
    
    sampling_freq = 1/100*(initi_freq*final_freq/numbr_freq)*ref_freq
    
    x = np.arange(0.0,10e-6,1/(spec_size*sampling_freq))[...,np.newaxis]
    x_size = np.size(x,axis=0)
    
    x_spec_full = np.arange(0.0,1/sampling_freq,1/(spec_size*sampling_freq))[...,np.newaxis]

    dom = np.linspace(-.5,.5,x_size)
    gauss = 1/(sigma * np.sqrt(2.0 * np.pi)) * np.exp( - (dom - mu)**2 / (2 * sigma**2) )[...,np.newaxis]

    central_range = ref_freq*(.48*dom + final_freq/(100*2))

    data_set,grid,respc,scale,ndx0 = load_dataset(conve_mod, gridname, output_path, filename)

    SP = Spectrum(initi_freq,final_freq,numbr_freq,parti_freq,x,gauss,spec_size)

    return SP, offset, dtsr, x_size, central_range, data_set, grid,respc,scale , ndx0 , grid_scale, x_spec_full

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

def set_domain_plot(grid,ndx0, data_set,grid_scale):
    
    this_series = np.empty(np.prod(grid.shape))
    this_series[ndx0] = np.nan
    this_series[np.logical_not(ndx0)] = np.zeros(data_set['doma'].shape[1])
    
    this_grid = grid*grid_scale

    fg,axs = plt.subplots(1,1,figsize=(10,7))
    h0 = axs.pcolormesh(this_grid.real,this_grid.imag,this_series.reshape(this_grid.shape),shading='nearest',vmin=-3,vmax=3)
    axs.set_aspect('equal', 'box')
    axs.tick_params(axis='x',which='both',bottom=True,top=False,labelbottom=True)
    axs.set_title('Total Field in Domain')

    fg.canvas.draw()
    fg.canvas.flush_events()

    return this_series, fg, h0

def create_matrix(SP, dtsr, x_size, central_range, data_set, ndx0, x_spec_full, mat_tseries, off, vlim, rlim, vlimmax, this_series, fg, h0):
    
    for i in range(0+1*int(1*x_size*4/8),1+1*int(1*x_size*4/8),1):
        _,_,spec = SP.synth_fseries_from_centr_freq(central_range[i-1])
        spec0 = spec[0:int(SP.spec_size)]
        
        for j in range(data_set['resp'].shape[-1]):
            resp_data = data_set['resp'][...,j]
            
            for k in range(data_set['resp'].shape[-2]):
                resptt = resp_data[...,k].transpose()
                new_resp = SP.expand_resp_new(resptt)
                tseries = SP.synth_tseries_from_spec_full_new(new_resp*spec0)            
                vec_tseries = tseries.transpose()
                vlimmax[k,j] = np.max(np.abs(vec_tseries),axis=-1)

                if k==j:
                    vlim[j,...] = vec_tseries
                if k==(j+off):
                    rlim[j,...] = vec_tseries                    


            if dtsr:
                doma_data = data_set['doma'][...,j]
                            
                for k in range(int(doma_data.shape[-1])):
                    domatt = doma_data[...,k].transpose()
                    new_doma = SP.expand_resp(domatt)
                    tseries = SP.synth_tseries_from_spec_full_new(new_doma*spec0)
                    mat_tseries[...,k] = tseries.transpose()

                mlim = np.amax(np.abs(mat_tseries))/4

                for ll in range(int(SP.spec_size*2/dtsr)):
                    this_series[np.logical_not(ndx0)] = mat_tseries[(ll+1)*dtsr-1,:].real/mlim

                    h0.set_array(this_series.ravel())
                    
                    fg.suptitle('time-step ' + str((ll+1)*dtsr-1) + '| time ' + str(x_spec_full[(ll+1)*dtsr-1]) + ' | height-step ' + str(j))
                    fg.canvas.draw()
                    fg.canvas.flush_events()

def plot_save_table(x_spec_full, off, vlim, rlim, vlimmax, probv):
    
    spec_grid,resp_grid = np.meshgrid(x_spec_full,probv)
    tx_grid,rx_grid = np.meshgrid(probv,probv)

    x_grid = (tx_grid*np.cos(np.pi/4) - rx_grid*np.sin(np.pi/4))/np.sqrt(2)
    y_grid = (tx_grid*np.sin(np.pi/4) + rx_grid*np.cos(np.pi/4))/np.sqrt(2)

    fg,axs = plt.subplots(2,1,figsize=(10,7))
 
    axs[0].set_title('Diagonal Signal')
    axs[0].pcolormesh(spec_grid,resp_grid,np.abs(vlim),shading='nearest')

    vlimmaxN = vlimmax/np.max(vlimmax)

    t_grid = np.diagonal(y_grid,offset=off)
    n_grid = np.diagonal(x_grid,offset=off)
    tpec_grid,tesp_grid = np.meshgrid(x_spec_full,t_grid)
    
    axs[1].set_title('{}-Off Diagonal Signal'.format(probv[int(probv.size/2)+off]))
    axs[1].pcolormesh(tpec_grid,tesp_grid,np.abs(rlim),shading='nearest')

    plt.figure(figsize=(7,7))
    plt.pcolormesh(x_grid,y_grid,vlimmaxN,shading='nearest',cmap='Greys')

    x_table = probv
    y_table = np.diagonal(vlimmaxN)
    y_off_table = np.diagonal(vlimmaxN,offset=off)

    d_table = np.array([x_table,y_table]).transpose()
    d_off_table = np.array([t_grid,y_off_table]).transpose()

    plt.figure('off_diagonal_matrix')
    plt.plot(x_table,y_table)
    plt.plot(t_grid,y_off_table)
    plt.title('Off Diagonal Matrix')
    plt.tight_layout()

    np.savetxt('vlimmax.csv',d_table, delimiter=',', header=','.join(('t','s')), comments='')
    np.savetxt('vlimmax_off.csv',d_off_table, delimiter=',', header=','.join(('t','s')), comments='')

def set_empty_matrix(SP, offset, data_set, respc, scale):
    mat_tseries = np.empty((SP.spec_size,data_set['doma'].shape[-2]),dtype='complex')
        
    vlim = np.empty((data_set['resp'].shape[-1],SP.spec_size*2),dtype='complex')
    rlim = np.empty(((data_set['resp'].shape[-1]-offset),SP.spec_size*2),dtype='complex')
    
    vlimmax = np.empty((int(data_set['resp'].shape[-2]),int(data_set['resp'].shape[-1])),dtype='float')

    probv =  respc+(.5+np.arange(-int(data_set['resp'].shape[-1])/2,int(data_set['resp'].shape[-1])/2))[:]*scale
    return mat_tseries,vlim,rlim,vlimmax,probv

def tseries(pre_config, set_empty_matrix, set_domain_plot, create_matrix, plot_save_table, config_file, output_path):

    SP, offset, dtsr, x_size, central_range, data_set, grid,respc,scale , ndx0, grid_scale, x_spec_full = pre_config(config_file,output_path)

    mat_tseries, vlim, rlim, vlimmax, probv = set_empty_matrix(SP, offset, data_set, respc, scale)

    this_series, fg, h0 = set_domain_plot(grid,ndx0,data_set,grid_scale)

    create_matrix(SP, dtsr, x_size, central_range, data_set, ndx0, x_spec_full, mat_tseries, offset, vlim, rlim, vlimmax, this_series, fg, h0)

    plot_save_table(x_spec_full, offset, vlim, rlim, vlimmax, probv)

    plt.show()


if __name__ == '__main__':
   
    if len(sys.argv) != 3:
        raise ValueError('Invalid number of arguments. Usage: {} output_path config.yaml'.format(sys.argv[0]))

    output_path = sys.argv[1]
    config_file = sys.argv[2]

    tseries(pre_config, set_empty_matrix, set_domain_plot, create_matrix, plot_save_table, config_file, output_path)