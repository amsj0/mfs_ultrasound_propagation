from pathlib import Path
import numpy as np
import h5py

def load_para(pathname,modifier,gridname):
    
    f = h5py.File(pathname + 'P' + modifier + gridname + '.h5','r')

    gridx = f['D']['X']

    gridy = f['D']['Z']

    shape = f['D']['s']

    radis = f['nRD']

    lambd = f['g']['wav']

    elemw = f['g']['eleme_wav']

    intec = f['g']['interf_centre']
   
    pistc = f['g']['piston_centre']

    grid = np.empty(gridx.shape,dtype='complex')
    
    grid.real = gridx;grid.imag = gridy;grid.shape = shape

    scale = (lambd[()],(elemw[()]/100))

    respc = (pistc[()]- intec[()])*(scale[0]/scale[1])
    
    grid =  grid*lambd[()]/radis[()]

    return respc,scale


def load_(pathname,para,filename):
    
    f =  h5py.File(pathname + para + '_' +  filename,'r')
        
    return f[para]

def save_dict_to_hdf5(M, path, datafile):
    with h5py.File(path + datafile + '.h5','w') as f:
        for item, dict in M.items():
            try:
                f_dict = dict.__dict__
                f_item = f.create_group(item)
                for k, v in f_dict.items():
                    f_item[k] = v   
            except Exception:
                f[item] = dict
    

def save_data_to_hdf5(data, path, attrs):

    # Assume: d_time_table, m_time_table, d_table, vlimmaxN, d_off_table already exist
    # and conve_mod, freq, output_path are defined.

    conve_mod = attrs['conve_mod']
    freq = attrs['freq']
    kr = 10*attrs.get('kr', None)
    dr = 10*attrs.get('dr', None))
    scale = attrs.get('wavelength', None)
    
    
    out_dir = Path(path)
    out_dir.mkdir(parents=True, exist_ok=True)
    h5_path = out_dir / f"bundle_{conve_mod}_{freq}_{kr}_{dr}.h5"

    with h5py.File(h5_path, "w") as h5:
        # File-level metadata (optional but handy)
        h5.attrs["conve_mod"] = str(conve_mod)
        h5.attrs["freq"] = str(freq)
        h5.attrs["kr"] = str(kr)
        h5.attrs["dr"] = str(dr)
        h5.attrs["wavelength"] = str(scale)

        def write_dataset(name, arr, dimensions=None, vectors=None):
            arr = np.asarray(arr)
            dset = h5.create_dataset(
                name,
                data=arr,
                compression="gzip",        # compact
                compression_opts=4,
                shuffle=True,              # better compression
                fletcher32=True            # integrity check
            )
            if dimensions is not None:
                # store column names like your CSV headers
                dset.attrs["dimensions"] = np.array(dimensions, dtype="S")
            if vectors is not None:
                for i, ax in enumerate(vectors):
                    dset.attrs[f"axis_{i}"] = ax

        for item, entry in data.items():
            arr    = entry.arr
            header = entry.header
            axis   = entry.axis
            write_dataset(item, arr, dimensions=header, vectors=axis)

        # Example usage:

    print(f"Saved: {h5_path}")


def create_keysized_to_hdf5(key, size, path, dataset):
    with h5py.File(path + dataset + '.h5', 'w') as f:
        f.create_dataset(key, size, dtype='complex')

def save_keyvalue_to_hdf5(key, value, path, dataset):
    with h5py.File(path + dataset + '.h5', 'w') as f:
        f[key] = value

def append_keyvalue_to_hdf5(key, value, ndx, path, dataset):
    with h5py.File(path + dataset + '.h5', 'a') as f:
        f[key][ndx,...] = value