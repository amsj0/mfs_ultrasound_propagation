import numpy as np
import h5py

class Store:
    
    def __init__(self, datafile, dataroot, metadata = None):
        
        self.Solution = { # Solution
            datafile : None
        }
        self.Behavior = { # Behavior
            datafile : None
        }
        self.Configuration = { # Configuration
            dataroot : metadata
        }

    def load_dict_to_store_solution(self, M, datafile):
        
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
                self.Solution[datafile][item] = np.copy(dict)

    def load_dict_to_store_behaviour(self, M, ndx):
        
        for key, value in M.items():
            self.Behavior[key][ndx,...] = value
        print('DataFile {} appended'.format(key))
       
    def init_store_behavior(self, datafile, size):
        
        if '' in self.Behavior:
            self.Behavior.pop('')

        self.Behavior[datafile] = np.empty(size,dtype='complex')