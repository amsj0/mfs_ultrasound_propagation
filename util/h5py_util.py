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
    kr = int(10*attrs.get('kr', None))
    dr = int(10*attrs.get('dr', None))
    scale = attrs.get('wavelength', None)
    
    
    out_dir = Path(path)
    out_dir.mkdir(parents=True, exist_ok=True)
    h5_path = out_dir / f"data_{conve_mod}_{freq}_{kr}_{dr}.h5"

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

# ----------------- Toeplitz helpers -----------------
def _toeplitz_symmetric_from_lastrow(last_row: np.ndarray) -> np.ndarray:
    """
    Rebuild an m×m symmetric Toeplitz matrix T from its last row r.
    Steps:
      1) v = r[::-1]  gives lag vector [t0, t1, ..., t_{m-1}]
      2) s = [t_{m-1}, ..., t1] + [t0, t1, ..., t_{m-1}]  (palindrome extension, len=2m-1)
      3) T[i,j] = s[i - j + (m - 1)]
    """
    r = np.asarray(last_row)
    if r.ndim != 1:
        raise ValueError("last_row must be 1D.")
    m = r.shape[0]
    v = r[::-1]
    s = np.concatenate([v[:0:-1], v])  # [v_{m-1}..v1] + [v0..v_{m-1}]  -> len 2m-1
    idx = np.arange(m)[:, None] - np.arange(m)[None, :] + (m - 1)
    return s[idx]


def _extract_block_lastrows(MT: np.ndarray):
    """
    Split MT into four m×m blocks and return their last rows.
    Blocks: TL=A, TR=B, BL=C, BR=D.
    """
    if MT.ndim != 2 or MT.shape[0] != MT.shape[1]:
        raise ValueError("MT must be a square 2D matrix.")
    N = MT.shape[0]
    if N % 2 != 0:
        raise ValueError("N must be even so the matrix can be split into 4 equal blocks.")
    m = N // 2

    A = MT[:m, :m]
    B = MT[:m, m:]
    C = MT[m:, :m]
    D = MT[m:, m:]

    return N, m, A[-1, :].copy(), B[-1, :].copy(), C[-1, :].copy(), D[-1, :].copy()


def _rebuild_from_block_lastrows(N: int, rA: np.ndarray, rB: np.ndarray, rC: np.ndarray, rD: np.ndarray, dtype=None):
    """
    Rebuild full N×N with symmetric Toeplitz sub-blocks from their last rows.
      A = toeplitz_symmetric_from_lastrow(rA), etc.
    """
    if N % 2 != 0:
        raise ValueError("N must be even to rebuild 4 equal blocks.")
    m = N // 2
    for name, r in (("rA", rA), ("rB", rB), ("rC", rC), ("rD", rD)):
        if r.ndim != 1 or r.shape[0] != m:
            raise ValueError(f"{name} must be a 1D vector of length N/2.")

    A = _toeplitz_symmetric_from_lastrow(rA)
    B = _toeplitz_symmetric_from_lastrow(rB)
    C = _toeplitz_symmetric_from_lastrow(rC)
    D = _toeplitz_symmetric_from_lastrow(rD)

    out_dtype = dtype or np.result_type(A, B, C, D)
    MT = np.empty((N, N), dtype=out_dtype)
    MT[:m, :m] = A.astype(out_dtype, copy=False)
    MT[:m, m:] = B.astype(out_dtype, copy=False)
    MT[m:, :m] = C.astype(out_dtype, copy=False)
    MT[m:, m:] = D.astype(out_dtype, copy=False)
    return MT


# ---- (Optional) compatibility with prior "diagonals" format ----
def _toeplitz_symmetric_from_diagvec(diagvec: np.ndarray) -> np.ndarray:
    v = np.asarray(diagvec)
    if v.ndim != 1:
        raise ValueError("diagvec must be 1D.")
    m = v.shape[0]
    idx = np.arange(m)
    return v[np.abs(idx[:, None] - idx[None, :])]

def _rebuild_from_block_diagvecs(N: int, dA: np.ndarray, dB: np.ndarray, dC: np.ndarray, dD: np.ndarray, dtype=None):
    A = _toeplitz_symmetric_from_diagvec(dA)
    B = _toeplitz_symmetric_from_diagvec(dB)
    C = _toeplitz_symmetric_from_diagvec(dC)
    D = _toeplitz_symmetric_from_diagvec(dD)
    out_dtype = dtype or np.result_type(A, B, C, D)
    m = N // 2
    MT = np.empty((N, N), dtype=out_dtype)
    MT[:m, :m] = A.astype(out_dtype, copy=False)
    MT[:m, m:] = B.astype(out_dtype, copy=False)
    MT[m:, :m] = C.astype(out_dtype, copy=False)
    MT[m:, m:] = D.astype(out_dtype, copy=False)
    return MT


# ----------------- I/O API (drop-in) -----------------
def load_hdf5_to_array(datafile, path, key):
    """
    Loads:
      - NEW: group with format='block_toeplitz_lastrow' (rA,rB,rC,rD) and rebuilds via palindrome extension.
      - LEGACY: dataset at f[key] (plain array), or group with format='block_toeplitz_diagonals' (dA..dD).
    """
    with h5py.File(path + datafile + '.h5', 'r') as f:
        obj = f[key]
        # Legacy plain dataset
        if isinstance(obj, h5py.Dataset):
            return np.copy(obj[...])

        if isinstance(obj, h5py.Group):
            fmt = obj.attrs.get('format', '')
            if isinstance(fmt, bytes):
                fmt = fmt.decode()

            if fmt == 'block_toeplitz_lastrow':
                N = int(obj.attrs['N'])
                dtype_str = obj.attrs.get('dtype', None)
                dtype = np.dtype(dtype_str) if dtype_str is not None else None

                rA = obj['rA'][...]
                rB = obj['rB'][...]
                rC = obj['rC'][...]
                rD = obj['rD'][...]
                return _rebuild_from_block_lastrows(N, rA, rB, rC, rD, dtype=dtype)

            # Back-compat with the previous diagonals format
            if fmt == 'block_toeplitz_diagonals':
                N = int(obj.attrs['N'])
                dtype_str = obj.attrs.get('dtype', None)
                dtype = np.dtype(dtype_str) if dtype_str is not None else None

                dA = obj['dA'][...]
                dB = obj['dB'][...]
                dC = obj['dC'][...]
                dD = obj['dD'][...]
                return _rebuild_from_block_diagvecs(N, dA, dB, dC, dD, dtype=dtype)

        raise ValueError(f"Unsupported HDF5 structure at key '{key}'.")


def save_dict_to_hdf5_toeplitz(M, path, datafile, *, diag_type: str = "anti", store_from_transpose: bool = False):
    """
    Saves a dict to HDF5.
    For 'MT' (square array), store one diagonal per sub-block:
      diag_type='anti' -> anti-diagonals; rebuilt as symmetric Hankel blocks.
      diag_type='main' -> main diagonals; rebuilt as symmetric Toeplitz blocks.
    If store_from_transpose=True, we extract diagonals from MT.T and mark that in metadata.
    """
    with h5py.File(path + datafile + '.h5', 'w') as f:
        for item, obj in M.items():
            if item == 'MT' and isinstance(obj, np.ndarray) and obj.ndim == 2 and obj.shape[0] == obj.shape[1]:
                N, m, rA, rB, rC, rD = _extract_block_lastrows(obj)
                g = f.create_group(item)
                g.attrs['format'] = 'block_toeplitz_lastrow'
                g.attrs['N'] = N
                g.attrs['dtype'] = str(obj.dtype)

                for name, vec in (('rA', rA), ('rB', rB), ('rC', rC), ('rD', rD)):
                    g.create_dataset(name, data=vec, compression='gzip', compression_opts=4, shuffle=True)
                continue

            # Generic fallback (legacy-like)
            try:
                f_obj = getattr(obj, '__dict__', None)
                f_item = f.create_group(item)
                if isinstance(f_obj, dict):
                    for k, v in f_obj.items():
                        try:
                            f_item.create_dataset(k, data=v)
                        except TypeError:
                            f_item.attrs[k] = v
                else:
                    if isinstance(obj, np.ndarray):
                        f.create_dataset(item, data=obj, compression='gzip', compression_opts=4, shuffle=True)
                    else:
                        ds = f.create_group(item)
                        ds.attrs['value'] = obj
            except Exception:
                if isinstance(obj, np.ndarray):
                    f.create_dataset(item, data=obj, compression='gzip', compression_opts=4, shuffle=True)
                else:
                    g = f.create_group(item)
                    g.attrs['value'] = str(obj)

    
    
def create_keysized_to_hdf5(key, size, path, dataset):
    with h5py.File(path + dataset + '.h5', 'w') as f:
        f.create_dataset(key, size, dtype='complex')

def save_keyvalue_to_hdf5(key, value, path, dataset):
    with h5py.File(path + dataset + '.h5', 'w') as f:
        f[key] = value

def append_keyvalue_to_hdf5(key, value, ndx, path, dataset):
    with h5py.File(path + dataset + '.h5', 'a') as f:
        f[key][ndx,...] = value