import numpy as np

vcos = np.frompyfunc(np.cos,1,1)
vsin = np.frompyfunc(np.sin,1,1)


class structtype():
    pass

def fn_copy_filter(seed, *copies): 

    for copy in copies:
        
        mask = np.logical_not(copy.ndx)
        
        copy.x = seed.x[mask]
        copy.y = seed.y[mask]
        copy.z = seed.z[mask]
        copy.a = seed.a[mask]
        copy.p = seed.p[mask]

def fn_enclosure_rectan(coor = None,sign = None,thres = None,disc = None): 

    varargout = []
    for c,s,t in zip(coor,sign,thres):
        varargout.append(s * c <= s * (t + s * disc ))
    
    return varargout

def fn_discretize_geometry_domain(G,centre_vector,grid_ratio):
    
    M = structtype()
    
    vec = []
    siz = []
    
    
    for grid,centre in zip(G,centre_vector):
        
        grid_number = int(np.floor(grid / grid_ratio))      
        grid_centre = grid_number - 0.5

        siz.append(grid_number * 2 + 1 * 0)
        vec.append(centre + (np.arange(-grid_centre,grid_centre+1)) * grid_ratio)
        
    M.x = np.repeat(vec[0][:,np.newaxis],siz[0],axis=1).reshape(-1)
    M.z = np.repeat(vec[1][np.newaxis,:],siz[1],axis=0).reshape(-1)
    
    M.X = M.x
    M.Z = M.z

    M.y = np.zeros((np.prod(siz)))
    M.a = np.zeros((np.prod(siz)))
    M.p = np.zeros((np.prod(siz)))
    M.s = siz

    return M

def fn_discretize_geometry_plane(RN = None,Centre = None, Orientation = None,Nelt = None,mode = None):

    P = structtype() 
    
    N = 1

    radi = 0
    pts_number = int(np.floor(RN / 2 * Nelt) * 2)
    RN = (pts_number - 1) / Nelt
        
    if pts_number > 1:
        radi += np.linspace(- RN / 2,RN / 2,pts_number)
    
    thet =              np.zeros((pts_number,))
    area = 1 / Nelt   * np.ones((pts_number,))
    norm = - Orientation * np.ones((pts_number,))

    #curv = np.inf * np.ones((1,pts_number))
    #assignin('base','cm',np.amin(curva))
    
    M = pts_number * N
    radiC = []
    areaC = []
    normC = []
    thetC = []

    for k in np.arange(N).reshape(-1):
        radiC.append( radi )
        areaC.append( area )
        normC.append( norm + k * (2 * np.pi / N) )
        thetC.append( thet + k * (2 * np.pi / N) )
    
    radi = np.array(radiC , dtype=object)
    thet = np.array(thetC , dtype=object)
    area = np.array(areaC , dtype=object)
    norm = np.array(normC , dtype=object)
    
    P.a = area
    P.x = (radi * vsin(thet) + Centre[0] )
    P.z = (radi * vcos(thet) + Centre[1] )
    P.y = P.z * 0
    P.n = norm

    P.a = np.hstack(P.a).astype(float)
    P.x = np.hstack(P.x).astype(float)
    P.z = np.hstack(P.z).astype(float)
    P.y = np.hstack(P.y).astype(float)
    P.n = np.hstack(P.n).astype(float)

    return P

def fn_surface_rectangular_stacking(P,gap = 0.5):
    
    # RECTANGULAR SHIFT BASED ON RECTANGULAR PACKING

    f0 = 1.0
    
    f = []
    f.append( f0 - gap)
    f.append( f0 + gap)
      
    c = f0 * np.exp(-1j*P.n)
    co = c * (f[1] - f0)
    ci = c * (f0 - f[0])
    
    P.co = P.c + co
    P.ci = P.c - ci

def f_surface_circular_stacking(P,Np,lambda0,gap = 0.0015):

    f0 = 1.0

    s = 3
    a = np.ceil(np.log((1 - np.sqrt(3) * np.pi/ Np) ** s - gap * lambda0) / np.log((1 - np.sqrt(3) * np.pi/ Np)))
    hap = ((1 + np.sqrt(3) * np.pi/ Np) ** a - 1)
    a = - np.ceil(np.log((1 - np.sqrt(3) * np.pi/ Np) ** s - gap * lambda0) / np.log((1 - np.sqrt(3) * np.pi/ Np)))
    ham = - (1 - (1 - np.sqrt(3) * np.pi/ Np) ** a)
    
    #ff = (f0 - 2 * np.sin(np.pi/ Np) / (1 + np.sin(np.pi/ Np)))
    #f = []
    #f.append(ff ** (s))
    #f.append(ff ** (- s))
   
    P.co = P.c * (f0 - ham)
    P.ci = P.c * (f0 + hap)