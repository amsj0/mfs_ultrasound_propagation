#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
import h5py
import time
import sys


# In[ ]:


spec_size = int(1/2*512)
hspc_size = int(spec_size/2)
#central_freq = 19/(.6*2) #15.833
#central_freq = 15.707963268 #5pi
central_freq = .5e6
range_freq = .1e6
ref_freq = 1e6

numbr_freq = 100
initi_freq = 1
final_freq = 100
parti_freq = 100
converge ='A'
modifier = 'r'

display_time_step_ratio = 16
dtsr = display_time_step_ratio

filename = '_' + str(numbr_freq) + '_' + str(initi_freq) + '_' + str(final_freq) + '.h5'
# pathname = '../octave/'
pathname = sys.argv[1]

print(pathname)

mu, sigma = 0, 0.1


# In[ ]:

sampling_freq = 1/100*(initi_freq*final_freq/numbr_freq)*ref_freq

print(sampling_freq)


factor = 1

x = np.arange(0.0,10e-6,1/(spec_size*sampling_freq))[:,None]
x_size = np.size(x,axis=0)
x_spec = np.arange(0.0,1/sampling_freq,1/(spec_size*sampling_freq))[:,None]
x_spec_full = np.arange(0.0,factor/sampling_freq,1/(spec_size*sampling_freq))[:,None]
dom = np.linspace(-.5,.5,x_size)
gauss = 1/(sigma * np.sqrt(2.0 * np.pi)) * np.exp( - (dom - mu)**2 / (2 * sigma**2) )[:,None]
#gauss = 1
central_range = range_freq*dom + central_freq


# In[ ]:


def synth_fseries_from_centr_freq(cent_freq):
    osc = np.exp(1j*cent_freq*2*np.pi*x)
    spec = np.fft.fft(osc*gauss,n=spec_size,axis=0)
    freq = np.fft.fftfreq(spec_size)

    return osc,freq,spec


# In[ ]:


def load_para(modifier):
    
    f = h5py.File(pathname + 'PP' + modifier + filename,'r')

    gridx = f['b']['value']['A']['value']['x']['value']

    gridy = f['b']['value']['A']['value']['z']['value']

    grid = np.empty(gridx.shape,dtype='complex')

    grid.real = gridx;grid.imag = gridy
    
    return grid/f['nRD']['value']


# In[ ]:


def load_resp(converge):
    
    f = h5py.File(pathname + 'resp_enh_' + str(converge) +  filename,'r')

    resp = np.array(f['resp']['value']).view(complex)

    return resp


# In[ ]:


def load_doma(converge):
    
    f = h5py.File(pathname + 'doma_enh_' + str(converge) +  filename,'r')

    doma = np.array(f['doma']['value']).view(complex)

    return doma


# In[ ]:


def load_(para,converge):
    
    f = h5py.File(pathname + para + '_enh_' + str(converge) +  filename,'r')

    dataset = f[para]['value']
    
    list_keys = list(dataset.keys())
    
    dims = tuple(dataset[list_keys[-1]])

    return dataset,dims


# In[ ]:


def expand_resp(resp):
    leading_size = int(initi_freq/(final_freq-initi_freq)*(numbr_freq-1))
    leading_zero = np.zeros((leading_size,1))
    trailing_size = int(hspc_size-parti_freq-leading_size)
    trailing_zero = np.zeros((trailing_size,1))
    new_resp = np.empty((hspc_size,1),dtype='complex')
    flip_resp = np.flip(resp,axis=0)
    new_resp = np.vstack((
        leading_zero,
        resp[:,np.newaxis],
        trailing_zero)
    )

    return new_resp


# In[ ]:


def expand_resp_0(resp):
    leading_size = int(initi_freq/(final_freq-initi_freq)*(numbr_freq-1))
    leading_zero = np.zeros((leading_size,1))
    trailing_size = int(hspc_size-parti_freq-leading_size)
    trailing_zero = np.zeros((trailing_size,1))
    new_resp = np.empty((spec_size,1),dtype='complex')
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


# In[ ]:


def expand_resp_1(resp):
    leading_size = int(initi_freq/(final_freq-initi_freq)*(numbr_freq-1))
    leading_zero = np.zeros((leading_size,resp.shape[1]))
    #print(leading_size)
    trailing_size = int(hspc_size-parti_freq-leading_size)
    trailing_zero = np.zeros((trailing_size,resp.shape[1]))
    #print(trailing_size)
    new_resp = np.empty((spec_size,resp.shape[1]),dtype='complex')
    flip_resp = np.flip(resp,axis=0)
    #print(2*(trailing_size+leading_size+resp.shape[0]))
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


# In[ ]:


def expand_resp_2(resp):
    leading_size = int(initi_freq/(final_freq-initi_freq)*(numbr_freq-1))
    leading_zero = np.zeros((leading_size,resp.shape[1],resp.shape[2]))

    trailing_size = int(hspc_size-leading_size-parti_freq)
    trailing_zero = np.zeros((trailing_size,resp.shape[1],resp.shape[2]))

    new_resp = np.empty((spec_size,resp.shape[1],resp.shape[2]),dtype='complex')

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


# In[ ]:


def synth_tseries(new_resp,cent_freq):
    osc = np.sin(cent_freq*np.pi*x)

    spec = np.fft.fft(osc*gauss,n=spec_size,axis=0)
    freq = np.fft.fftfreq(spec_size)

    new_full = new_resp*spec[:, np.newaxis]
    tseries = np.fft.ifft(new_full,n=spec_size,axis=0)
    return osc,new_full,tseries


# In[ ]:


def synth_tseries_from_spec(new_resp,spec):
    new_full = new_resp*spec
    tseries = np.fft.ifft(new_full,n=spec_size,axis=0)
    return new_full,tseries


# In[ ]:


def synth_tseries_from_spec_full(new_resp,spec,factor):
    new_full = new_resp*spec
    tseries = np.fft.ifft(new_full,n=int(factor*spec_size),axis=0)
    return new_full,tseries


# In[ ]:


def plt_arrange(ii):
    
    osc,new_full,tseries = synth_tseries(central_range[ii])
    #axs = plt.figure().subplots(3,1)
    '''
    ax = fig.add_subplot(gs[0,0])
    ax.plot(x,osc*gauss)
    ax.set_title('pulse')

    ax = fig.add_subplot(gs[:,2])
    ax.imshow(abs(new_full[:,0,:]))
    ax.set_title('spectrum vs distance')

    ax = fig.add_subplot(gs[:,1])
    ax.imshow(tseries[:,0,:].real)
    ax.set_title('a-scan vs distance')
    
    ax = fig.add_subplot(gs[:,0])
    ax.imshow(abs(new_resp))
    ax.set_title('spectrum vs distance')
    
    ax = fig.add_subplot(gs[1,0])
    ax.plot(np.max(abs(tseries),axis=0))
    ax.set_title('max vs distance')
    '''
    plt.plot(np.max(abs(tseries),axis=0).transpose())
    #plt.show()
    plt.savefig('tseries' + str(ii) + '.png')
    plt.clear(True)


# In[ ]:


def plt_max_tseries(converge):
    resp = load_resp(converge)
    new_resp = expand_resp_2(resp)
    max_tseries = np.empty((x_size,new_resp.shape[2],new_resp.shape[1]),dtype='double')
    #mat_tseries = np.empty((x_size,new_resp.shape[0],new_resp.shape[1],new_resp.shape[2]),dtype='complex')
     
    for ii in range(x_size):
        osc,new_full,tseries = synth_tseries(new_resp,central_range[ii])
        max_tseries[ii,:,:] = np.max(abs(tseries),axis=0).transpose()
        #mat_tseries[ii,:,:,:] = tseries
    #return mat_tseries
    return max_tseries


# In[ ]:


def plt_mat_tseries_1(converge):
    resp = load_resp(converge)
    
    #max_tseries = np.empty((x_size,new_resp.shape[2],new_resp.shape[1]),dtype='double')
    #mat_tseries = np.empty((x_size,spec_size,resp.shape[1],resp.shape[2]),dtype='complex')
    mat_tseries = np.empty((x_size,resp.shape[2]),dtype='complex')
     
    for ii in range(int(x_size)):
        osc,freq,spec = synth_fseries_from_centr_freq(central_range[ii])
        for jj  in range(int(resp.shape[2])):
            new_resp = expand_resp_1(resp[:,:,jj])
            new_full,tseries = synth_tseries_from_spec(new_resp,spec)
            tsumseries = np.sum(tseries,axis=-1)
            tstdseries = np.std(tsumseries)
        #max_tseries[ii,:,:] = np.max(abs(tseries),axis=0).transpose()
            mat_tseries[ii,jj] = tstdseries
    return mat_tseries
    #return max_tseries


# In[ ]:


#get_ipython().run_line_magic('matplotlib', 'notebook')


# In[ ]:


doma_dataset,doma_dims = load_('doma',converge + modifier)

doma_list_keys = list(doma_dataset.keys())
doma_data = np.array(doma_dataset[doma_list_keys[0]]['value']).view(complex)

resp_dataset,resp_dims = load_('resp',converge + modifier)

resp_list_keys = list(resp_dataset.keys())
resp_data = np.array(resp_dataset[resp_list_keys[0]]['value']).view(complex)

grid = load_para(modifier)


# In[ ]:


# Testing t-domain analysis on grid
this_time_ndx = 0

plt.ion()

#fg = plt.figure(figsize=(10,7))
fg,axs = plt.subplots(4,1,figsize=(10,7))

#this_series = 20*np.log10(abs(mat_tseries[this_time_ndx,:,2]))
this_series = np.zeros((1,doma_data.shape[1]))
#this_series = mat_tseries[this_time_ndx,:,2].real
#this_series = abs(mat_tseries[this_time_ndx,:,2])
#this_series = this_doma[:,0]
this_grid = grid
this_series.shape = this_grid.shape
h0 = axs[0].pcolormesh(this_grid.real,this_grid.imag,this_series,shading='nearest',vmin=-1, vmax=1)
h1 = axs[1].pcolormesh(this_grid.real,this_grid.imag,this_series,shading='nearest',vmin=-1, vmax=1)
h2 = axs[2].pcolormesh(this_grid.real,this_grid.imag,this_series,shading='nearest',vmin=-1, vmax=1)
h3 = axs[3].pcolormesh(this_grid.real,this_grid.imag,this_series,shading='nearest',vmin=-1, vmax=1)
axs[0].set_aspect('equal', 'box')
axs[0].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
axs[0].set_title('Refracted')
axs[1].set_aspect('equal', 'box')
axs[1].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
axs[1].set_title('Reflected')
axs[2].set_aspect('equal', 'box')
axs[2].tick_params(axis='x',which='both',bottom=False,top=False,labelbottom=False)
axs[2].set_title('Incident')
axs[3].set_aspect('equal', 'box')
axs[3].set_title('Total')
#plt.colorbar(h0,cax=axs[0])
fg.canvas.draw()
fg.canvas.flush_events()

print(int(doma_dims[0]))

mat_tseries = np.empty((int(factor*spec_size),doma_data.shape[1],4),dtype='complex')
#factor = 1.25
for ii in range(0+1*int(1*x_size*4/8),1+1*int(1*x_size*4/8),1):
    osc,freq,spec = synth_fseries_from_centr_freq(central_range[ii-1])
    spec0 = spec[0:hspc_size]
    for jj in range(0,int(doma_dims[0]),int((doma_dims[0]-1)/5)):
        index = np.ravel_multi_index((jj,0),doma_dims)
        doma_data = np.array(doma_dataset[doma_list_keys[index]]['value']).view(complex)
        index = np.ravel_multi_index((jj,1),doma_dims)
        doma_data0 = np.array(doma_dataset[doma_list_keys[index]]['value']).view(complex)
        index = np.ravel_multi_index((jj,2),doma_dims)
        doma_data1 = np.array(doma_dataset[doma_list_keys[index]]['value']).view(complex)
        doma_data2 = doma_data0 + doma_data1 + doma_data
        #domat = doma[:,:,jj*8]
        for kk in range(int(doma_data.shape[1])):
            domatt = doma_data[:,kk].transpose()
            new_doma = expand_resp(domatt)
            new_full,tseries = synth_tseries_from_spec_full(new_doma,spec0,factor)
            mat_tseries[:,kk,0] = tseries.transpose()
            domatt = doma_data0[:,kk].transpose()
            new_doma = expand_resp(domatt)
            new_full,tseries = synth_tseries_from_spec_full(new_doma,spec0,factor)            
            mat_tseries[:,kk,1] = tseries.transpose()
            domatt = doma_data1[:,kk].transpose()
            new_doma = expand_resp(domatt)
            new_full,tseries = synth_tseries_from_spec_full(new_doma,spec0,factor)            
            mat_tseries[:,kk,2] = tseries.transpose()
            domatt = doma_data2[:,kk].transpose()
            new_doma = expand_resp(domatt)
            new_full,tseries = synth_tseries_from_spec_full(new_doma,spec0,factor)            
            mat_tseries[:,kk,3] = tseries.transpose()             
        vlim = np.amax(np.abs(mat_tseries))/4
        for ll in range(int(factor*int(spec_size/dtsr))):
            this_series = mat_tseries[(ll+1)*dtsr-1,:,0].real/vlim
            #this_series = 20*np.log10(abs(mat_tseries[ll*dtsr,:,0]))
            this_series.shape = this_grid.shape
            h0.set_array(this_series.ravel())
            
            this_series = mat_tseries[(ll+1)*dtsr-1,:,1].real/vlim
            #this_series = 20*np.log10(abs(mat_tseries[ll*dtsr,:,1]))
            this_series.shape = this_grid.shape            
            h1.set_array(this_series.ravel())
            
            this_series = mat_tseries[(ll+1)*dtsr-1,:,2].real/vlim
            #this_series = 20*np.log10(abs(mat_tseries[ll*dtsr,:,2]))
            this_series.shape = this_grid.shape            
            h2.set_array(this_series.ravel())
            
            this_series = mat_tseries[(ll+1)*dtsr-1,:,3].real/vlim
            #this_series = 20*np.log10(abs(mat_tseries[ll*dtsr,:,2]))
            this_series.shape = this_grid.shape            
            h3.set_array(this_series.ravel())
            
            
            fg.suptitle('time-step ' + str((ll+1)*dtsr-1) + '| time ' + str(x_spec_full[(ll+1)*dtsr-1]) + ' | height-step ' + str(jj))
            #domat.shape

            #plt.draw()
            #time.sleep(1e-4)
            #clear_output(wait=False)
            #plt.show()#draw();
            fg.canvas.draw()
            #plt.grid(True
            fg.canvas.flush_events()
