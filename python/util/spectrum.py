import numpy as np

class Spectrum:
        
    def __init__(self,initi,final,parti,numbr,x,gauss,spec_size):

        self.initi = initi
        self.final = final
        self.parti = parti
        self.numbr = numbr

        self.x = x
        self.gauss = gauss
        self.spec_size = int(spec_size/2)


    def expand_resp_new(self,resp):
        leading_size = int(self.initi/(self.final-self.initi)*(self.numbr-1))
        trailing_size = int(self.spec_size-self.parti-leading_size)
        new_resp = np.zeros((self.spec_size,1),dtype='complex')
        new_resp[leading_size:-trailing_size] = resp[:,np.newaxis]
        
        return new_resp

    def expand_resp(self,resp):
        leading_size = int(self.initi/(self.final-self.initi)*(self.numbr-1))
        leading_zero = np.zeros((leading_size,1))
        trailing_size = int(self.spec_size-self.parti-leading_size)
        trailing_zero = np.zeros((trailing_size,1))
        new_resp = np.empty((self.spec_size,1),dtype='complex')
        flip_resp = np.flip(resp,axis=0)
        new_resp = np.vstack((
            leading_zero,
            resp[:,np.newaxis],
            trailing_zero)
        )

        return new_resp

    def expand_resp_0(self,resp):
        leading_size = int(self.initi/(self.final-self.initi)*(self.numbr-1))
        leading_zero = np.zeros((leading_size,1))
        trailing_size = int(self.spec_size-self.parti-leading_size)
        trailing_zero = np.zeros((trailing_size,1))
        new_resp = np.empty((self.spec_size,1),dtype='complex')
        flip_resp = np.flip(resp,axis=0)
        new_resp.real = np.vstack((
            leading_zero,
            resp[:,np.newaxis].real,
            trailing_zero,
            trailing_zero,
            flip_resp[:,np.newaxis].real,
            leading_zero)
        )

        new_resp.imag = np.vstack((
            leading_zero,
            resp[:,np.newaxis].imag,
            trailing_zero,
            trailing_zero,
            -flip_resp[:,np.newaxis].imag,
            leading_zero)
        )
        return new_resp


    def expand_resp_1(self,resp):
        leading_size = int(self.initi/(self.final-self.initi)*(self.numbr-1))
        leading_zero = np.zeros((leading_size,resp.shape[1]))
        trailing_size = int(self.spec_size-self.parti-leading_size)
        trailing_zero = np.zeros((trailing_size,resp.shape[1]))
        new_resp = np.empty((self.spec_size,resp.shape[1]),dtype='complex')
        flip_resp = np.flip(resp,axis=0)
        new_resp.real = np.vstack((
            leading_zero,
            resp.real,
            trailing_zero,
            trailing_zero,
            flip_resp.real,
            leading_zero)
        )

        new_resp.imag = np.vstack((
            leading_zero,
            resp.imag,
            trailing_zero,
            trailing_zero,
            -flip_resp.imag,
            leading_zero)
        )
        return new_resp


    def expand_resp_2(self,resp):
        leading_size = int(self.initi/(self.final-self.initi)*(self.numbr-1))
        leading_zero = np.zeros((leading_size,resp.shape[1],resp.shape[2]))

        trailing_size = int(self.spec_size-leading_size-self.parti)
        trailing_zero = np.zeros((trailing_size,resp.shape[1],resp.shape[2]))

        new_resp = np.empty((self.spec_size,resp.shape[1],resp.shape[2]),dtype='complex')

        new_resp.real = np.vstack((
            leading_zero,
            resp.real,
            trailing_zero,
            trailing_zero,
            resp[::-1].real,
            leading_zero)
        )

        new_resp.imag = np.vstack((
            leading_zero,
            resp.imag,
            trailing_zero,
            trailing_zero,
            -resp[::-1].imag,
            leading_zero)
        )
        return new_resp


    def synth_fseries_from_centr_freq(self,cent_freq):
        osc = np.exp(1j*2*np.pi*(cent_freq*self.x+1/4))*self.gauss
        spec = np.fft.fft(osc,n=self.spec_size*2,axis=0)
        freq = np.fft.fftfreq(self.spec_size*2)

        return osc,freq,spec

    def synth_tseries(self,new_resp,cent_freq):
        osc = np.sin(cent_freq*np.pi*self.x)
        spec = np.fft.fft(osc*self.gauss,n=self.spec_size*2,axis=0)
        freq = np.fft.fftfreq(self.spec_size*2)

        new_full = new_resp*spec[:, np.newaxis]
        tseries = np.fft.ifft(new_full,n=self.spec_size*2,axis=0)
        return osc,new_full,tseries

    def synth_tseries_from_spec(self,new_resp,spec):
        new_full = new_resp*spec
        tseries = np.fft.ifft(new_full,n=self.spec_size*2,axis=0)
        return new_full,tseries


    def synth_tseries_from_spec_full(self,new_resp,spec):
        new_full = new_resp*spec
        tseries = np.fft.ifft(new_full,n=self.spec_size*2,axis=0)
        return new_full,tseries


    def synth_tseries_from_spec_full_new(self,spec):
        tseries = np.fft.ifft(spec,n=self.spec_size*2,axis=0)
        return tseries