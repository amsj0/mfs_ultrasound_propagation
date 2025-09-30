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

    return np.block([[sol_ST, sol_TT]])

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