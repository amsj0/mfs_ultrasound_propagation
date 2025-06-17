import numpy as np

def heuristic(nRD, g):
    
    # WAVENUMBER VECTOR RANGE
    k0 = 2*np.pi*np.linspace(nRD,nRD,g.nff)*(g.sfr/g.fr)
    kr = g.skr/1000
    k0_length = k0.shape[0]
    kr_length = kr.shape[0]
    # DENSITY VECTOR RANGE
    d0 = np.array(g.rj)
    dr = g.sdr/1000
    d0_length = d0.shape[0]
    dr_length = dr.shape[0]
    # SENSITIVITY VECTOR RANGE
    p0 = np.ones(1)
    pr = g.spr/100
    p0_length = p0.shape[0]
    pr_length = pr.shape[0]
    # WAVENUMBER VECTOR RANGE

    sfr = g.sfr/1000

    RD = nRD*g.wav
    max_k0 = np.amax(k0) # Maximum wavenumber

    # SET DIMENSIONS FOR RECTANGULAR GRID
    lambda0 = 1
    d_cur = d0
    ac_iomega = 2*np.pi*g.sfr
    rho_r = lambda rd: g.rj[0]*rd
    spd_r = lambda rc: g.cj[0]*rc
    delta_diff = lambda rd: (4*g.att/3)/rho_r(rd)
    c_c = lambda rc,rd: (spd_r(rc)**2-1j*ac_iomega*delta_diff(rd))**0.5
    rho_var = lambda rc,rd: spd_r(rc)**2/(c_c(rc,rd))**2
    g.cjR = lambda rc,rd: c_c(rc,rd)
    g.rjR = lambda rc,rd: (rho_r(rd)*rho_var(rc,rd))
    g.keq = lambda rc,rd: np.sqrt(-(1j*ac_iomega/g.cjR(rc,rd))**2)

    return k0,kr,kr_length,dr,dr_length,sfr,RD,lambda0,d_cur