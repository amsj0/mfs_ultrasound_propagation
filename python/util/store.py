import numpy as np
import h5py

class Store:
    
    def __init__(self, datafile):
        
        self.Solution = { # Solution
            datafile : None
        }

    def load_dict_to_hdf5(self, M, datafile):
        
        if '' in self.Solution:
            self.Solution.pop('')

        self.Solution[datafile] = {}
        #f = h5py.File(datafile,'w', driver='core', backing_store=False, block_size = 12800000)
        for item, dict in M.items():
            try:
                f_dict = dict.__dict__
                f_item = self.Solution[datafile][item]
                for k, v in f_dict.items():
                    f_item[k] = v   
            except Exception:
                self.Solution[datafile][item] = dict