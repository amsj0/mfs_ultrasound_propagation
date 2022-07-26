import numpy as np
import h5py

def load_para(pathname,modifier,gridname):
    
    f = h5py.File(pathname + 'P' + modifier + gridname + '.h5','r')

    gridx = f['D']['X']

    gridy = f['D']['Z']

    ndx0 = f['D']['ndx0']

    shape = f['D']['s']

    radis = f['nRD']

    lambd = f['g']['wav']

    elemw = f['g']['eleme_wav']

    intec = f['g']['interf_centre']
   
    pistc = f['g']['piston_centre']

    grid = np.empty(gridx.shape,dtype='complex')
    
    grid.real = gridx;grid.imag = gridy;grid.shape = shape

    scale = lambd[()]/(elemw[()]/100)

    respc = (pistc[()]- intec[()])*scale
    
    grid =  grid*lambd[()]/radis[()]

    return grid,respc,scale,ndx0


def load_(pathname,para,filename):
    
    f =  h5py.File(pathname + para + '_' +  filename,'r')
        
    return f[para]

def save_dict_to_hdf5(M, datafile):
    with h5py.File(datafile + '.h5','w') as f:
        for item, dict in M.items():
            try:
                f_dict = dict.__dict__
                f_item = f.create_group(item)
                for k, v in f_dict.items():
                    f_item[k] = v   
            except Exception:
                f[item] = dict

def create_keysized_to_hdf5(key, size, path, dataset):
    with h5py.File(path + dataset + '.h5', 'w') as f:
        f.create_dataset(key, size, dtype='complex')

def save_keyvalue_to_hdf5(key, value, path, dataset):
    with h5py.File(path + dataset + '.h5', 'w') as f:
        f[key] = value

def append_keyvalue_to_hdf5(key, value, ndx, path, dataset):
    with h5py.File(path + dataset + '.h5', 'a') as f:
        f[key][ndx,...] = value