import numpy as np
from numpy.linalg import inv
from scipy.linalg import solve

def inc_ref(p1,v1,p2,v2):

    v1 = inv(v1)
    p1 = p1 @ v1
    
    v1 = inv(p2-p1 @ v2)

    p2_invS = p2 @ v1
    v2_invS = v2 @ v1

    p1 = inv(v2-inv(p1) @ p2)

    p2_invT = p2 @ p1
    v2_invT = v2 @ p1

    #return np.concatenate((np.concatenate((p2 @ v1,p2 @ p1),axis=1),np.concatenate((v2 @ v1,v2 @ p1),axis=1)),axis=0)
    
    #return np.block([[b2_invS,b2_invT]])
    return np.block([[p2_invS, p2_invT],
                    [v2_invS, v2_invT]])

def inc_ref2(p1, v1, p2, v2):
    """
    All inputs are square matrices with compatible shapes.
    This version avoids explicit inverses (np.linalg.inv) in favor of solves.
    """
    # Step 1: p1 := p1 @ inv(v1)  -> right-side solve: solve(v1.T, p1.T).T
    p1_r = np.linalg.solve(v1.T, p1.T).T
    # Step 2: S := p2 - p1 @ v2   (note: @ has higher precedence than -)
    S = p2 - p1_r @ v2
    T = v2 - np.linalg.solve(p1_r, p2)

    # We need p2 @ inv(S) and v2 @ inv(S) -> left solves on S^T
    RHS = np.concatenate((p2.T, v2.T), axis=1)  # shape: n x (2n)
    sol_ST = np.linalg.solve(S.T, RHS).T
    sol_TT = np.linalg.solve(T.T, RHS).T
    #p2_invS, v2_invS = np.split(sol_ST, 2, axis=0)
    #p2_invS = np.linalg.solve(S.T, p2.T).T
    #v2_invS = np.linalg.solve(S.T, v2.T).T
    # Step 3: T := v2 - inv(p1) @ p2  -> inv(p1) @ p2 = solve(p1, p2)
    # We need p2 @ inv(T) and v2 @ inv(T) -> left solves on T^T
    #p2_invT, v2_invT = np.split(sol_TT, 2, axis=0)
    #p2_invT = np.linalg.solve(T.T, p2.T).T
    #v2_invT = np.linalg.solve(T.T, v2.T).T
    # Assemble the 2x2 block:
    # [[p2 @ inv(S), p2 @ inv(T)],
    #  [v2 @ inv(S), v2 @ inv(T)]]
    return np.block([[sol_ST, sol_TT]])
    #return np.block([[p2_invS, p2_invT],
    #                 [v2_invS, v2_invT]])


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