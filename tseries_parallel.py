import os
os.environ.setdefault("OMP_NUM_THREADS", "1")
os.environ.setdefault("OPENBLAS_NUM_THREADS", "1")
os.environ.setdefault("MKL_NUM_THREADS", "1")
os.environ.setdefault("VECLIB_MAXIMUM_THREADS", "1")
os.environ.setdefault("NUMEXPR_NUM_THREADS", "1")

import numpy as np
from scipy.signal import hilbert
import scipy.signal as sp

from util.spectrum import Spectrum
from util.store import Store
from util.entry import Entry
from util.h5py_util import *
import yaml

from itertools import product
from concurrent.futures import ProcessPoolExecutor

import sys

from util.h5py_util import *

def load_dataset(task_name,store ,conve_mod, gridname, pathname, skr, sdr):
    
    filename = gridname + '_' + str(int(skr)) + '_' + str(int(sdr))

    respc,scale = load_para(pathname, conve_mod, gridname)
    
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

    return data_set,respc,scale,(skr,sdr)


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
            - scale (float): Scale factors from the loaded dataset.
            - ndx0 (int): Index value from the loaded dataset.
            - x_spec_full (numpy.ndarray): Full spectrum x array.
            - conve_mod (str): Conversion model from the configuration.
            - tranx_num (int): Number of transducers from the configuration.
    Note:
        freq_scale is calculated as the ratio of the specified spectrum size to the number of frequencies divided by the spectrum band.
    """
   
    if task_name == "tseries":
        with open(config_file, 'r') as f:
            cfg = yaml.safe_load(f)
            initi_freq = cfg['initi_freq']
            final_freq = cfg['final_freq']
            parti_freq = cfg['parti_freq']
            numbr_freq = cfg['numbr_freq']
            conve_mod = cfg['CONVE_MOD']

            sdr = [cfg['sdr']]
            skr = [cfg['skr']]
    else:

        dataroot = list(store.Configuration)[-1]

        parts = dataroot.split('_')
        try :
            final_freq = int(parts[-1])
            initi_freq = int(parts[-2])
            numbr_freq = int(parts[-3])
            conve_mod  = '_'.join(parts[:-3])

            parti_freq = numbr_freq

            skr = store.Configuration[dataroot].skr
            sdr = store.Configuration[dataroot].sdr
        except:
            print('Error in the name of the dataset stored in the Store object. The name must be in the format CONVE_MOD_NUMBR_FREQ_INITI_FREQ_FINAL_FREQ')
        
    with open(config_file, 'r') as f:
        cfg = yaml.safe_load(f)
        offset = cfg['offset']

        tranx_num = cfg['tranx_num']
        spec_band = cfg['spec_band']
        spec_size = cfg['spec_size']
        ref_freq = cfg['ref_freq']
        input_mod = cfg['INPUT_MOD']
            
              
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
    f_dom = np.linspace(1/x_size,1,x_size) - 0.5

    central_range = ref_freq*(1*f_dom + final_freq/(100*2))

    input = np.empty((x_size,1),dtype='complex')
    if input_mod == 'experiment':
        input = np.loadtxt(output_path + 'experimental_pulse.txt', dtype=np.complex_)[np.newaxis,:]
    elif input_mod == 'synthetic_gauss':
        input = np.sin(2 * np.pi * central_range * x_dom) * gauss_filter
    elif input_mod == 'synthetic_cosine':
        input = np.exp(1j * 2 * np.pi * (central_range * x - 1 / 5)) * cosine_filter.T
        #np.sin(2 * np.pi * central_range * x_dom) * filter
    else:
        raise ValueError("Invalid INPUT_MOD. Choose 'experiment' or 'synthetic' or 'synthetic_cosine'.")
    
    gridname = '_' + str(numbr_freq) + '_' + str(initi_freq) + '_' + str(final_freq) 

    SP = Spectrum(initi_freq,final_freq,numbr_freq,parti_freq,x,filter,input,spec_size)
    
    results = [load_dataset(task_name, store, conve_mod, gridname, output_path, kr, dr) for kr, dr in product(skr, sdr)]
    
    data_set,respc,scale,par_map = map(list, zip(*results))

    return SP, offset, x_size, central_range, data_set,respc,scale,par_map , x_spec_full, conve_mod, tranx_num


def _worker_j(args):
    (j, resp_j, spec0_col, SP,
     matrix_x_strafe, matrix_y_strafe, offset) = args

    data_y_size = resp_j.shape[-1]
    vlimmax_col = np.empty((data_y_size,), dtype=np.float64)
    vlim_j = None
    rlim_pair = None

    mlim_updates = {}  # (mx,my) -> vec_ts (last one wins within this j)

    matrix_x_ndx = j % matrix_x_strafe
    on_x_stride = (matrix_x_ndx == 0)

    for k in range(data_y_size):

        matrix_y_ndx = k % matrix_y_strafe
        on_y_stride = (matrix_y_ndx == 0)
        # your original extraction
        resptt = resp_j[..., k].T

        new_resp = np.conj(SP.expand_resp_new(resptt))
        tseries  = SP.synth_tseries_from_spec_full_new(new_resp * spec0_col)
        vec_ts   = tseries.T  # shape (1, 2*spec_size)

        # vlimmax column entry
        vlimmax_col[k] = np.max(np.abs(vec_ts), axis=-1)

        # vlim row when k == j
        if k == j:
            vlim_j = vec_ts

        # rlim row when k == j + offset
        if k == j + offset:
            ridx = j - int((abs(offset) - offset) / 2)
            rlim_pair = (ridx, vec_ts)

        # mlim sparse updates on stride hits
        if on_x_stride and on_y_stride:
            mx = int(-0.5 + j / matrix_x_strafe)
            my = int(-0.5 + k / matrix_y_strafe)
            mlim_updates[(mx, my)] = vec_ts  # keep last, matches serial overwrite

    # compress updates for return
    mlim_list = [(mx, my, arr) for (mx, my), arr in sorted(mlim_updates.items(), key=lambda x: x[0][1])]
    return j, vlim_j, vlimmax_col, rlim_pair, mlim_list

def run_per_j_processes(lim,vlimmax,SP,
                        data_set, spec0,
                        matrix_x_strafe, matrix_y_strafe, offset,
                        max_workers=None, chunksize=4):
    resp = np.asarray(data_set["resp"])  # expected shape (..., data_y_size, data_x_size)
    den = resp[..., 0, 0][..., np.newaxis, np.newaxis]
    np.divide(resp, den, out=resp, where=den != 0)
    
    #data_y_size = resp.shape[-2]
    data_x_size = resp.shape[-1]

    vlim,mlim,rlim = lim

    # build per-j tasks (send only one resp slice per task)
    tasks = []
    for j in range(data_x_size):
        resp_j = np.asarray(resp[..., j])            # shape (..., data_y_size)
        tasks.append((j, resp_j, spec0, SP,
                      matrix_x_strafe, matrix_y_strafe, offset))

    if max_workers is None:
        max_workers = os.cpu_count() or 4

    # process in parallel; map preserves input order (j ascending)
    with ProcessPoolExecutor(max_workers=max_workers) as ex:
        for j, vlim_j, vlimmax_col, rlim_pair, mlim_list in ex.map(_worker_j, tasks, chunksize=chunksize):
            # assemble results (in j order to match original semantics)
            vlimmax[:, j] = vlimmax_col
            if vlim_j is not None:
                vlim[j, ...] = vlim_j
            if rlim_pair is not None:
                ridx, arr = rlim_pair
                if 0 <= ridx < rlim.shape[0]:
                    rlim[ridx, ...] = arr
            # apply mlim updates for this j in ascending my order (same as serial k loop)
            for mx, my, arr in mlim_list:
                if 0 <= mx < mlim.shape[0] and 0 <= my < mlim.shape[1]:
                    mlim[mx, my, ...] = arr

    #return dict(vlim=vlim, rlim=rlim, mlim=mlim, vlimmax=vlimmax)


def create_matrix(SP, x_size, central_range, data_set, x_spec_full, mat_tseries, offset, lim, vlimmax):
    respttm = []

    
    vlim,mlim,rlim = lim


    data_x_size = data_set['resp'].shape[-1]
    data_y_size = data_set['resp'].shape[-2]

    matrix_x_strafe = (data_x_size-1)/(mlim.shape[0])
    matrix_y_strafe = (data_y_size-1)/(mlim.shape[1])

    i = x_size // 2 - 1
    central_freq = central_range[i]

    #osc,freq,spec = SP.synth_fseries_from_centr_freq(central_freq)
    osc,freq,spec = SP.synth_fseries_from_input(i)
    spec0 = spec[0:int(SP.spec_size)]
    

    print('Central frequency: ' + str(central_freq) + ' Hz' )

    run_per_j_processes(
        lim, vlimmax, SP,
        data_set,spec0,
        matrix_x_strafe=int(matrix_x_strafe), matrix_y_strafe=int(matrix_y_strafe), offset=offset,
        max_workers=(os.cpu_count()), chunksize=16
    )
        

    return int(central_freq/1000)
    

def save_table(conve_mod, x_spec_full, freq, par, lim, vlimmax, prob, offset, scale, output_path=''):
    
    vlimmaxN = vlimmax/np.max(vlimmax[0]) # Normalize the maximum value of the diagonal signal
    #vlimmaxN = vlimmax/np.max(vlimmax[-1])

    vlim,mlim,rlim = lim
    probv,probm = prob

    tx_grid,rx_grid = np.meshgrid(probv,probv)

    y_grid = (tx_grid*np.sin(np.pi/4) + rx_grid*np.cos(np.pi/4))/np.sqrt(2)

    t_grid = np.diagonal(y_grid,offset=offset)
    y_off_table = np.diagonal(vlimmaxN,offset=offset)
    
    d_table = np.diagonal(vlimmaxN)
    d_off_table = y_off_table
    d_time_table = np.real(vlim)
    m_time_table = np.real(mlim)

    data = {
        "M": Entry(vlimmaxN, ("wavelenght", "wavelenght"),(probv,probv)),
        "N": Entry(m_time_table, ("wavelenght","wavelenght", "time"),(probm[...,0],probm[...,1],x_spec_full.transpose())),
        "T": Entry(d_time_table, ("wavelenght", "time"),(probv,x_spec_full.transpose())),
        "V": Entry(d_table, ("wavelenght",),(probv,)),
        "Voff": Entry(d_off_table, ("wavelenght",),(t_grid,))
    }
    
    save_data_to_hdf5(
        data,
        output_path,
        {'conve_mod': conve_mod,'freq': freq, 'kr': par[0], 'dr': par[1], 'wavelength': scale[0]}
        )

def set_empty_matrix(SP, offset, doma_size, resp_size, respc, scale, tranx_num=1):
    

    mat_tseries = np.empty((SP.spec_size,doma_size[-2]),dtype='complex')
        
    data_m_size = resp_size[-2],resp_size[-1]
    respv_size = (np.array(data_m_size) - 1) / 2

    probmn = tuple()
    probvn = tuple()
    for respv_s in respv_size:
        respvn = 0 * respc + (np.arange(-respv_s, 1 + respv_s))[:] / scale[1]
        step_n = int(2 * respv_s / tranx_num)
        start_n = len(respvn) // 2 - (tranx_num // 2) * step_n + int(step_n / 2)
        stop_n = start_n + tranx_num * step_n
        probmn += (respvn[start_n:stop_n:step_n],)
        probvn += (respvn,)


    probt = np.meshgrid(*probmn)
    probm = np.stack(probt,axis=-1)
    probv = probvn[0]

    vlim = np.empty((probv.shape[0],SP.spec_size*2),dtype='complex')
    mlim = np.empty((probm.shape[0],probm.shape[1],SP.spec_size*2),dtype='complex')
    rlim = np.empty(((resp_size[-1]-np.abs(offset)),SP.spec_size*2),dtype='complex')
    
    vlimmax = np.empty((int(resp_size[-2]),int(resp_size[-1])),dtype='float')
    
    return mat_tseries,(vlim,mlim,rlim),vlimmax,(probv,probm)

def tseries_parallel(store,config_file, output_path,task_name):

    SP, offset, x_size, central_range, data_set ,respc ,scale,par_map , x_spec_full, conve_mod, tranx_num = pre_config(task_name, store, config_file,output_path)

    for set, s_respc, s_scale,par in zip(data_set, respc, scale,par_map):
        print('Loaded dataset with doma shape ' + str(set["doma"].shape) + ' and resp shape ' + str(set["resp"].shape))

        mat_tseries, lim, vlimmax, prob = set_empty_matrix(SP, offset, set["doma"].shape, set["resp"].shape, s_respc, s_scale, tranx_num)

        freq = create_matrix(SP, x_size, central_range, set, x_spec_full, mat_tseries, offset, lim, vlimmax)

        save_table(conve_mod, x_spec_full, freq, par, lim, vlimmax, prob, offset, scale, output_path)

if __name__ == '__main__':
   
    if len(sys.argv) != 3:
        raise ValueError('Invalid number of arguments. Usage: {} /path/to/output  config.yaml'.format(sys.argv[0]))

    output_path = sys.argv[1]
    config_file = sys.argv[2]

    store = Store('','')
    tseries_parallel(
        store=store,
        config_file=config_file,
        output_path=output_path,
        task_name ="tseries"
    )