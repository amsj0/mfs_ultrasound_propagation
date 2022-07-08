import numpy as np
from scipy.io import loadmat

def fn_analyze(path = None,datafile = None,configfile = None,analysisfile = None):
    loadmat(np.array([path,datafile]))
    loadmat(np.array([path,configfile,'.h5']))
    loadmat(np.array([path,analysisfile]))
    #   DEFINES PISTON INDEXES
    ppt_per_surface = 1 + int(np.floor(g.piston_radius * (Neltoverlambda / 100)))


    ndx.T = fn_integrate_indexing(ppt_per_surface,len(b.Ti.x),len(b.To.x))
    ndx.R = fn_integrate_indexing(ppt_per_surface,len(b.Ri.x),len(b.Ro.x))
    response.pitch.size = (len(T.a) - ppt_per_surface + 1)
    response.catch.size = (len(R.a) - ppt_per_surface + 1)
    field.range.ptt = np.zeros((response.pitch.size,len(b.D.c)))
    strength.range.ptt = cell(response.pitch.size,response.catch.size)
    response.range.pid = np.zeros((response.pitch.size,response.catch.size))
    response.range.prr = np.zeros((response.pitch.size,response.catch.size))
    response.range.prl = np.zeros((response.pitch.size,response.catch.size))
    response.range.ptt = np.zeros((response.pitch.size,response.catch.size))
    response.range.ptx = np.zeros((response.pitch.size,1))
    response.range.prx = np.zeros((response.catch.size,1))
    response.pitch.A = np.zeros((1,response.pitch.size))
    response.catch.A = np.zeros((1,response.pitch.size))
    response.pitch.I = T.z(not b.Ti.ndx ) / nRD
    response.pitch.O = T.z(not b.To.ndx ) / nRD
    response.catch.I = R.z(not b.Ri.ndx ) / nRD
    response.catch.O = R.z(not b.Ro.ndx ) / nRD
    response.ndx.RO = find(not b.Ro.ndx )
    response.ndx.RI = find(not b.Ri.ndx )
    response.ndx.TO = find(not b.To.ndx )
    response.ndx.TI = find(not b.Ti.ndx )

    ### TODO FIX ME
    ndx.tst = reshape(np.arange(1,(response.pitch.size * response.catch.size)+1),response.pitch.size,response.catch.size)
    n1,n2 = ndgrid(np.arange(1,response.pitch.size+1),np.arange(1,response.catch.size+1))
    n0 = np.transpose((np.mod(n2 + n1 + 1,response.pitch.size) + 1))
    ndx.mat = np.transpose(reshape(ndx.tst((n0 - 0) + (n1 - 1) * response.catch.size),response.pitch.size,response.catch.size))

    for tt in np.arange(1,response.pitch.size+1).reshape(-1):
        ndxI = np.array([ndx.T.catI[:,ppt_per_surface - 1 + tt]])
        ndxO = np.array([ndx.T.catO[:,ppt_per_surface - 1 + tt]])
        response.pitch.A[tt] = mean(np.array([response.pitch.I(ndxI),response.pitch.O(ndxO)]))
        ndx.T.i[tt] = ndxI
        ndx.T.o[tt] = ndxO

    for rr in np.arange(1,response.catch.size+1).reshape(-1):
        ndxI = np.array([ndx.R.catI[:,ppt_per_surface - 1 + rr]])
        ndxO = np.array([ndx.R.catO[:,ppt_per_surface - 1 + rr]])
        response.catch.A[rr] = mean(np.array([response.catch.I(ndxI),response.catch.O(ndxO)]))
        ndx.R.i[rr] = ndxI
        ndx.R.o[rr] = ndxO

    for tt in np.arange(1,response.pitch.size+1).reshape(-1):
        #   EXTRACT TRANSMITTER PISTON INDEXES
        ndx.TI = ndx.T.i[tt]
        ndx.TO = ndx.T.o[tt]
        resTO = response.ndx.TO(ndx.TO)
        resTI = response.ndx.TI(ndx.TI)
        #   FILTER FIELD PARAMETERS
        b.D.prr[b.D.ndx0,:] = energy_ratio * np.sum(Mcmb[1,1](:,ndx.TO), 2-1)
        b.D.prr[not b.D.ndx0 ,:] = np.sum(Mcmb[2,1](:,ndx.TI), 2-1)
        b.D.prl[b.D.ndx0,:] = np.sum(Mcmb[3,1](:,ndx.TI), 2-1)
        b.D.prl[not b.D.ndx0 ,:] = energy_ratio * np.sum(Mcmb[4,1](:,ndx.TO), 2-1)
        b.D.pid[b.D.ndx0,:] = np.sum(Mcmb[5,1](:,ndx.TI), 2-1)
        b.D.pid[not b.D.ndx0 ,:] = energy_ratio * np.sum(Mcmb[6,1](:,ndx.TO), 2-1)
        b.D.ptt = b.D.pid + b.D.prl + b.D.prr
        field.range.prr[tt,:] = b.D.prr
        field.range.prl[tt,:] = b.D.prl
        field.range.pid[tt,:] = b.D.pid

        b.T.prr[0] = energy_ratio * np.sum(Mcmb[1,2](:,ndx.TO), 2-1)
        b.T.prr[2] = np.sum(Mcmb[2,2](:,ndx.TI), 2-1)
        b.T.prl[0] = np.sum(Mcmb[3,2](:,ndx.TI), 2-1)
        b.T.prl[2] = energy_ratio * np.sum(Mcmb[4,2](:,ndx.TO), 2-1)
        b.T.pid[0] = np.sum(Mcmb[5,2](:,ndx.TI), 2-1)
        b.T.pid[2] = energy_ratio * np.sum(Mcmb[6,2](:,ndx.TO), 2-1)
        b.T.ptt[0] = (b.T.pid[0] + b.T.prl[0] + b.T.prr[0])
        b.T.ptt[2] = (b.T.pid[2] + b.T.prl[2] + b.T.prr[2])
        #   EXTRACT RECEIVER PISTON INDEXES
        for rr in np.arange(1,response.catch.size+1).reshape(-1):

            ndx.RI = ndx.R.i[rr]
            ndx.RO = ndx.R.o[rr]

            resRO = response.ndx.RO(ndx.RO)
            resRI = response.ndx.RI(ndx.RI)
            #   EXTRACT PISTON CENTRE
	    #   UPDATE FIGURES
            #   ANALYSE RESPONSE RANGES
            response.range.pid[tt,rr] = (sum(b.T.pid[0](ndx.RI)) + sum(b.T.pid[2](ndx.RO)))
            response.range.prr[tt,rr] = (sum(b.T.prr[0](ndx.RI)) + sum(b.T.prr[2](ndx.RO)))
            response.range.prl[tt,rr] = (sum(b.T.prl[0](ndx.RI)) + sum(b.T.prl[2](ndx.RO)))
        #     if(rr==floor(response_size/2))
        response.range.ptx[tt] = (sum(response.range.prx))

    response.range.ptt = response.range.prl + response.range.prr + response.range.pid

    rd.rr = cat(3,response.range.pid,response.range.prl,response.range.prr)
    rd.dd = cat(3,field.range.prr,field.range.prl,field.range.pid)
    return rd

#    return rd
