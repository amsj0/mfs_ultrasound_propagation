import numpy as np
from numpy.linalg import inv
from scipy.linalg import solve

def inc_ref(p1,v1,p2,v2):

    v1 = inv(v1)
    p1 = p1 @ v1
    
    v1 = inv(p2-p1 @ v2)
    p1 = inv(v2-inv(p1) @ p2)

    return np.concatenate((np.concatenate((p2 @ v1,p2 @ p1),axis=1),np.concatenate((v2 @ v1,v2 @ p1),axis=1)),axis=0)

def inc_ref3(z1,y1,p2,v2):
 
    v1 = p2 - z1 @ v2
    p1 = v2 - y1 @ p2
        
    pppv = solve(v1.conj().T,p2.conj().T).conj().T
    pvpp = solve(p1.conj().T,p2.conj().T).conj().T
    vppv = solve(v1.conj().T,v2.conj().T).conj().T
    vvpp = solve(p1.conj().T,v2.conj().T).conj().T
    
    return np.concatenate((np.concatenate((pppv,pvpp),axis=1),np.concatenate((vppv,vvpp),axis=1)),axis=0)

def impedance(p,v):
    
    #z = p / v
    #y = v / p
    return solve(v.conj().T,p.conj().T).conj().T,solve(p.conj().T,v.conj().T).conj().T    

def transfer(TF,Th):
      
    return np.subtract(np.identity(TF.shape[0]),TF) @ Th

def transferOLD(Thi,Tho,Tri,Tro,TFI,TFO,TIh,TOh):
    
    MH = []
    MR = []

    #'''
    Thi_TFI = Thi @ TFI
    Tho_TFI = Tho @ TFI
    Thi_TFO = Thi @ TFO
    Tho_TFO = Tho @ TFO

    MH.append(Thi @ TOh - Thi_TFO @ TOh)
    #MH.append(Thi_TFO @ TOh)
    MH.append(Tho @ TIh - Tho_TFI @ TIh)
    #MH.append(Tho_TFI @ TIh)
    MH.append(Thi @ TIh - Thi_TFO @ TIh)
    #MH.append(Thi @ TIh - Thi_TFI @ TIh)
    MH.append(Tho @ TOh - Tho_TFI @ TOh)
    #MH.append(Tho @ TOh - Tho_TFO @ TOh)
    #MH.append(p0kI)
    #MH.append(p0kO)

    Tri_TFI = Tri @ TFI
    Tro_TFI = Tro @ TFI
    Tri_TFO = Tri @ TFO
    Tro_TFO = Tro @ TFO  

    MR.append(Tri @ TOh - Tri_TFO @ TOh)
    #MR.append(Tri_TFO @ TOh)
    MR.append(Tro @ TIh - Tro_TFI @ TIh)
    #MR.append(Tro_TFI @ TIh)
    MR.append(Tri @ TIh - Tri_TFO @ TIh)
    #MR.append( Tri_TFO @ TIh)
    #MR.append(Tri @ TIh - Tri_TFI @ TIh)
    MR.append(Tro @ TOh - Tro_TFI @ TOh)
    #MR.append( Tro_TFI @ TOh)
    #MR.append(Tro @ TOh - Tro_TFO @ TOh)
    #MR.append(p0mI)
    #MR.append(p0mO)
    #'''
    
    return MH,MR    