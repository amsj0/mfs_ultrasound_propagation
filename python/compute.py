#from attr import define
import numpy as np
import pyopencl as cl

from util.propagator import inc_ref,transfer

import os

MCHARGE = 0
'''
'''
os.environ['PYOPENCL_COMPILER_OUTPUT'] = '1'
os.environ['PYOPENCL_self.ctx'] = '0'

class Compute:
        
    def __init__(self,P,T,S):
        
        self.Probe = P
        self.Tranx = T
        self.Surfc = S
        
        nT = np.size(self.Tranx['emitter'].c)
        nS = np.size(self.Surfc['collo'].c)
        nR = np.size(self.Probe['receiver'].c)
        nD = np.size(self.Probe['domain'].c)

        self.T = { # Transfer
            'domain'   : np.zeros((nD,2*nS), dtype=np.complex128), # to Domain from Surface
            'receiver' : np.zeros((nR,2*nS), dtype=np.complex128) # to Receiver from Surface
        }
        
        self.B = { # Boundary
            'emitter'  : np.zeros((2*nS,nT), dtype=np.complex128), # real incidence on collo Surface
            'upper' : np.zeros((nS,2*nS), dtype=np.complex128), # virtual incidence on upper Surface
            'lower' : np.zeros((nS,2*nS), dtype=np.complex128)  # virtual incidence on lower Surface
        }

        self.F = { # Incident Field
            'domain'   : np.zeros((nD,nT), dtype=np.complex128),  # on Domain
            'receiver' : np.zeros((nR,nT), dtype=np.complex128)   # on Receiver
        }

        self.M = { # Behaviour Field
            'domain'   : np.zeros((nD,nT), dtype=np.complex128),  # on Domain
            'receiver' : np.zeros((nR,nT), dtype=np.complex128)   # on Receiver
        }

        self.MUT = np.zeros((2*nS,nT), dtype=np.complex128) # Upper Transfer
        self.MLT = np.zeros((2*nS,nT), dtype=np.complex128) # Lower Transfer        

        self.MT = np.zeros((2*nS,2*nS), dtype=np.complex128) # Transfer Matrix

        self.nT = nT
        self.nR = nR
        self.nS = nS
        self.nD = nD

    def InitCL(self, DEVICE="CPU"):
       
        platforms = cl.get_platforms()

        self.ctx = cl.Context(
            dev_type=cl.device_type.ALL,
            properties=[(cl.context_properties.PLATFORM, platforms[0])])
        self.queue = cl.CommandQueue(self.ctx)
        self.mf    = cl.mem_flags

        with open('src/pyopencl-hankel-complex.cl', 'r') as code:
            self.pgr = cl.Program(self.ctx,code.read()).build("-Isrc/ -cl-strict-aliasing")

    def besselh(self,host,m):

        h0 = np.empty((host.size), dtype=np.complex128)
        h1 = np.empty((host.size), dtype=np.complex128)

        buff  = cl.Buffer(self.ctx, self.mf.READ_ONLY | self.mf.COPY_HOST_PTR, hostbuf=host.astype(np.complex128))

        h0_buff = cl.Buffer(self.ctx, self.mf.WRITE_ONLY, h0.nbytes)
        h1_buff = cl.Buffer(self.ctx, self.mf.WRITE_ONLY, h1.nbytes)
        
        completeEvent = self.pgr.hankel_01_complex(self.queue, host.shape, None, buff, h0_buff, h1_buff, np.int32(1))
        
        events = [completeEvent]
        
        cl.enqueue_copy(self.queue, h0, h0_buff, wait_for=events)
        cl.enqueue_copy(self.queue, h1, h1_buff, wait_for=events)
        
        h0.shape = host.shape
        h1.shape = host.shape
        
        bh0 = h0      
        bh1 = h1
        
        del h0_buff,h1_buff
        
        return bh0, bh1

    def propagator_side(self, R, side, C, p_cur, p_out ):

        m = MCHARGE
        kind = 0
        k_cur,d_cur = p_cur
        k_out,d_out = p_out        

        Ar = C.a[np.newaxis,:]

        factor = 1j/4
        
        mask = R.m[side]

        Ma = R.c[mask,np.newaxis] - C.c

        #midsize = Ma.shape[1]

        #Th = np.zeros((Ma.shape[0],2*Ma.shape[1]), dtype=np.complex128)

        nx = C.n[np.newaxis,:]

        rcos = (Ma.real*np.cos(nx)+Ma.imag*np.sin(nx))/np.abs(Ma)

        z_host = k_cur*np.abs(Ma)

        bh0, bh1 = self.besselh(z_host,m)   
        
        #Th[:,midsize:] = -1j/4*area[np.newaxis,:]*bh0
        #Th[:,:midsize] = 1j/4*area[np.newaxis,:]*bh1*rcos*k_cur/k_r*dr
        
        return  factor*Ar*bh1*rcos*(k_cur/d_cur),factor*Ar*bh0 # testing

    def field_boundary(self, R, C, side , p_cur, p_out):

        m = MCHARGE
        kind = 0
        k_cur,d_cur = p_cur
        k_out,d_out = p_out

        mask = C.m[side]

        Ar = C.a[np.newaxis,mask]

        Ma = R.c[:,np.newaxis] - C.c[mask]
        
        #midsize = Ma.shape[0]
        #Th = np.zeros((2*midsize,Ma.shape[1]), dtype=np.complex128)

        ny = R.n[:,np.newaxis]

        
        rcos = (Ma.real*np.cos(ny)+Ma.imag*np.sin(ny))/np.abs(Ma)
        rsin = (Ma.real*np.sin(ny)-Ma.imag*np.cos(ny))/np.abs(Ma)
        
        ex0 = np.exp(1j*(m*np.angle(Ma)-C.n[np.newaxis,mask]))
        
        z_host = k_cur*np.abs(Ma)
        
        bh0, bh1 = self.besselh(z_host,m)
            
        #Th[:midsize,:] = area[:,np.newaxis]*1*bh0*ex0*k_r
        #Th[midsize:,:] = area[:,np.newaxis]*(1*bh1*rcos*k_cur + 1j*bh0/np.abs(Ma)*rsin*m)*ex0
        
        return Ar*1*bh0*ex0*d_cur,Ar*(1*bh1*rcos*k_cur + 1j*bh0/np.abs(Ma)*rsin*m)*ex0 # testing


    def field_side_m(self, R, side, p_cur, p_out ):

        m = MCHARGE
        kind = 0
        k_cur,d_cur = p_cur
        k_out,d_out = p_out

        #len_R = R.c.shape[0]

        Ma = R['collo'].c - R[side].c[:,np.newaxis]

        ny = R['collo'].n.transpose()

        rcos1 = (Ma.real*np.cos(ny)+Ma.imag*np.sin(ny))/np.abs(Ma)

        z_host = k_cur*np.abs(Ma)
        
        bh0, bh1 = self.besselh(z_host,m)
        
        #Bh = bh0
        #Dh = bh1*rcos1*k_cur/k_r*dr

        return bh0,bh1*rcos1*(k_cur/d_cur) # testing


    def reference(self,R,C,side,p_cur):

        m = MCHARGE
        kind = 0
        k_cur = p_cur[0]
        
        maskR = R.m[side]

        maskC = C.m[side]

        Ar = C.a[np.newaxis,maskC]

        Ma = R.c[maskR,np.newaxis]-C.c[maskC]
        
        z_host = k_cur*np.abs(Ma)

        bh0, _ = self.besselh(z_host,m)

        return bh0*Ar

    def compute_field_upper_side_m(self,S,p_cur,p_out):

        self.B['upper'][...,:self.nS],self.B['upper'][...,self.nS:] = self.field_side_m(S,'upper',p_cur,p_out)

    def compute_field_lower_side_m(self,S,p_cur):

        self.B['lower'][...,:self.nS],self.B['lower'][...,self.nS:] = self.field_side_m(S,'lower',p_cur,p_cur)

    def compute_propagator_upper_side(self,P,probe,side,S,p_cur,p_out):

        self.T[probe][P.m[side],:self.nS],self.T[probe][P.m[side],self.nS:] = self.propagator_side(P,side,S['collo'],p_cur,p_out)

    def compute_propagator_lower_side(self,P,probe,side,S,p_cur):

        self.T[probe][P.m[side],:self.nS],self.T[probe][P.m[side],self.nS:] = self.propagator_side(P,side,S['collo'],p_cur,p_cur)

    def compute_field_upper_boundary(self,S,T,side,p_cur,p_out):
        
        self.B['emitter'][:self.nS,T.m[side]],self.B['emitter'][self.nS:,T.m[side]] = self.field_boundary(S['collo'],T,side,p_cur,p_out)
    
    def compute_field_lower_boundary(self,S,T,side,p_cur):
        
        self.B['emitter'][:self.nS,T.m[side]],self.B['emitter'][self.nS:,T.m[side]] = self.field_boundary(S['collo'],T,side,p_cur,p_cur)

    def compute_reference(self,P,probe,T,side,p_cur):
        
        #self.F[probe][side[probe][...,np.newaxis]*side['emitter']] = self.reference(P,side[probe],T,side['emitter'],k_cur).reshape(-1)
        self.F[probe][P.m[side][...,np.newaxis]*T.m[side]] = self.reference(P,T,side,p_cur).reshape(-1)

    def compute_lower_side(self,p_cur):
        
        side = 0

        self.compute_field_lower_side_m(self.Surfc,p_cur)
        self.compute_field_lower_boundary(self.Surfc, self.Tranx['emitter'], side, p_cur)
        for probe, P in self.Probe.items():
            self.compute_propagator_lower_side(P,probe,side,self.Surfc,p_cur)
            self.compute_reference(P,probe,self.Tranx['emitter'],side,p_cur)

    def compute_upper_side(self,p_cur,p_out):
        
        side = 1

        self.compute_field_upper_side_m(self.Surfc,p_cur,p_out)
        self.compute_field_upper_boundary(self.Surfc, self.Tranx['emitter'], side, p_cur,p_out)
        for probe, P in self.Probe.items():
            self.compute_propagator_upper_side(P,probe,side,self.Surfc,p_cur,p_out)
            self.compute_reference(P,probe,self.Tranx['emitter'],side,p_cur)

    def propagate_lower_incref(self):

        p1 = self.B['upper'][:,:self.nS]
        v1 = self.B['upper'][:,self.nS:]
        p2 = -self.B['lower'][:,:self.nS]
        v2 = -self.B['lower'][:,self.nS:]

        self.MT = inc_ref(p1,v1,p2,v2)

    def propagate_upper_incref(self):

        p1 = self.B['lower'][:,:self.nS]
        v1 = self.B['lower'][:,self.nS:]
        p2 = -self.B['upper'][:,:self.nS]
        v2 = -self.B['upper'][:,self.nS:]

        self.MT = inc_ref(p1,v1,p2,v2)

    def propagate_incref(self,side):

        p1 = self.B['lower'][:,:self.nS]
        v1 = self.B['lower'][:,self.nS:]
        p2 = -self.B['upper'][:,:self.nS]
        v2 = -self.B['upper'][:,self.nS:]

        self.MT = inc_ref(p1,v1,p2,v2)
        if ~side:
            self.MT = np.identity(2*self.nS) - self.MT

    def propagate_transfer(self):
        
        self.MUT = transfer(self.MT,self.B['emitter'])
        self.MLT = transfer(self.MT,self.B['emitter'])

    def propagate_scatter(self):
        
        emitter_mask = self.Tranx['emitter'].m[0]
        nmitter_mask = self.Tranx['emitter'].m[1]

        for probe, Transfer in self.T.items():
            
            trans_mask = self.Probe[probe].m[0]
            nrans_mask = self.Probe[probe].m[1]
            probe_mask = trans_mask[...,np.newaxis]
            nrobe_mask = nrans_mask[...,np.newaxis]

            self.M[probe][probe_mask*emitter_mask] = (Transfer[trans_mask,...] @ self.MUT[...,emitter_mask]).flatten()
            self.M[probe][nrobe_mask*nmitter_mask] = -(Transfer[nrans_mask,...] @ self.MLT[...,nmitter_mask]).flatten()
            self.M[probe][probe_mask*nmitter_mask] = -(Transfer[trans_mask,...] @ self.MUT[...,nmitter_mask]).flatten()
            self.M[probe][nrobe_mask*emitter_mask] = (Transfer[nrans_mask,...] @ self.MLT[...,emitter_mask]).flatten()
            
            self.M[probe] += self.F[probe]
            self.M[probe][(probe_mask|nrobe_mask)*(emitter_mask&nmitter_mask)] /= 2

        return self.M

    def null(self):
        return 0