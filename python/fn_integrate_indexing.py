import numpy as np
    
def fn_integrate_indexing(points = None,varargin = None): 
    size_I = varargin[0]
    size_O = varargin[2]
    # MATRIX MAGIC
    matrix_ndx = lambda v = None: cellfun(lambda x = None: circshift(v,x,2),num2cell(np.arange(0,(points - 1)+1)),'UniformOutput',0)
    vec = np.array([num2cell(np.arange(1,size_I+1),1),cell(1,size_O)])
    mat = matrix_ndx(vec)
    ndx.catI = cat(1,mat[:])
    vec = np.array([cell(1,size_I),num2cell(np.arange(1,size_O+1),1)])
    mat = matrix_ndx(vec)
    ndx.catO = cat(1,mat[:])
    return ndx