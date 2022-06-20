#!/usr/bin/env python
# coding: utf-8

# In[ ]:


#from nbformat import read
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.animation as animation
#from scipy.signal import hilbert
import yaml
import h5py
import time
import sys


# In[ ]:


spec_size = int(2/2*512)
hspc_size = int(spec_size/2)
#central_freq = 19/(.6*2) #15.833
#central_freq = 15.707963268 #5pi

yaml_path = 'config.yaml'

if len(sys.argv) != 1:
    yaml_path = sys.argv[1]

with open(yaml_path, 'r') as f:
    cfg = yaml.safe_load(f)
    numbr_freq = cfg['numbr_freq']
    initi_freq = cfg['initi_freq']
    final_freq = cfg['final_freq']
    parti_freq = cfg['parti_freq']
    converge = cfg['converge']
    modifier = cfg['modifier']
    skr = cfg['skr']
    sdr = cfg['sdr']
    piston_radius = cfg['piston_radius']
    pathname = cfg['pathname']



ref_freq = 1e6

central_freq = ref_freq/100*final_freq/2
range_freq = ref_freq*.48

display_time_step_ratio = 256
dtsr = display_time_step_ratio

gridname = '_' + str(numbr_freq) + '_' + str(initi_freq) + '_' + str(final_freq) 

filename = gridname + '_' + str(sdr) + '_' + str(skr) + '_' + str(piston_radius) + '.h5'
# pathname = '../octave/'
# pathname = sys.argv[1]

print(pathname)

mu, sigma = 0, 0.1


# In[ ]:

sampling_freq = 1/100*(initi_freq*final_freq/numbr_freq)*ref_freq

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
    
    f = h5py.File(pathname + 'P' + modifier + gridname + '.h5','r')

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

def expand_resp_new(resp):
    leading_size = int(initi_freq/(final_freq-initi_freq)*(numbr_freq-1))
    #leading_zero = np.zeros((leading_size,1))
    trailing_size = int(hspc_size-parti_freq-leading_size)
    #trailing_zero = np.zeros((trailing_size,1))
    new_resp = np.zeros((hspc_size,1),dtype='complex')
    new_resp[leading_size:-trailing_size] = resp[:,np.newaxis]
    
    return new_resp

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

def synth_tseries_from_spec_full_new(spec,factor):
    tseries = np.fft.ifft(spec,n=int(factor*spec_size),axis=0)
    return tseries


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

if __name__ == '__main__':


    doma_dataset,doma_dims = load_('doma',converge + modifier)

    doma_list_keys = list(doma_dataset.keys())
    doma_data = np.array(doma_dataset[doma_list_keys[0]]['value']).view(complex)

    resp_dataset,resp_dims = load_('resp',converge + modifier)

    resp_list_keys = list(resp_dataset.keys())
    resp_data = np.array(resp_dataset[resp_list_keys[0]]['value']).view(complex)

    grid = load_para(converge + modifier)

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

    #print(int(doma_dims[0]))
    #print(int(np.sqrt(resp_dims[0])))

    mat_tseries = np.empty((int(factor*spec_size),doma_data.shape[1],4),dtype='complex')
    vec_tseries = np.empty((int(factor*spec_size),4),dtype='complex')
    
    #vlim = np.empty((int(resp_dims[0]),int(factor*spec_size)),dtype='complex')
    vlim = np.empty((int(resp_dims[0]),int(factor*spec_size)),dtype='complex')
    rlim = np.empty((int(resp_dims[0]),int(factor*spec_size)),dtype='complex')

    #vlim = np.empty((int(resp_dims[0]),resp_data.shape[1]),dtype='complex')
    #vlim = np.empty((int(doma_dims[0]),resp_data.shape[1]),dtype='complex')
    
    #factor = 1.25
    fg,ax = plt.subplots(1,1,figsize=(10,7))

    for ii in range(0+1*int(1*x_size*4/8),1+1*int(1*x_size*4/8),1):
        osc,freq,spec = synth_fseries_from_centr_freq(central_range[ii-1])
        spec0 = spec[0:hspc_size]
        
        for jj in range(0,int(resp_dims[0])):
        #for jj in range(0,int(resp_dims[0]),1+int(np.sqrt(resp_dims[0]))):

            index = np.ravel_multi_index((jj,0),resp_dims)
            
            resp_index = resp_dataset[resp_list_keys[index]]['value']

            if np.issubdtype(np.float64,resp_index.dtype):
                resp_data0 = np.array(resp_index).view(float)
            else:
                resp_data0 = np.array(resp_index).view(complex)

            index = np.ravel_multi_index((jj,1),resp_dims)          
            resp_index = resp_dataset[resp_list_keys[index]]['value']

            if np.issubdtype(np.float64, resp_index.dtype):
                resp_data1 = np.array(resp_index).view(float)
            else:
                resp_data1 = np.array(resp_index).view(complex)

            index = np.ravel_multi_index((jj,2),resp_dims)           
            resp_index = resp_dataset[resp_list_keys[index]]['value']
            
            if np.issubdtype(np.float64, resp_index.dtype):
                resp_data2 = np.array(resp_index).view(float)
            else:
                resp_data2 = np.array(resp_index).view(complex)
            
            resp_dataR = resp_data2 + resp_data1
            resp_data3 = resp_dataR + resp_data0
            
            """
            resptt = resp_data[:,0]
                                 
            new_resp = expand_resp(resptt)
            new_full,tseries = synth_tseries_from_spec_full(new_resp,spec0,factor)            
            vec_tseries[:,0] = tseries.transpose()
            
            resptt = resp_data0[:,0]
            new_resp = expand_resp(resptt)
            new_full,tseries = synth_tseries_from_spec_full(new_resp,spec0,factor)            
            vec_tseries[:,1] = tseries.transpose()
            
            resptt = resp_data1[:,0]
            new_resp = expand_resp(resptt)
            new_full,tseries = synth_tseries_from_spec_full(new_resp,spec0,factor)            
            vec_tseries[:,2] = tseries.transpose()
            """

            for kk in [jj]:
                    
                resptt = resp_dataR[:,kk].transpose()
                new_resp = expand_resp_new(resptt)
                tseries = synth_tseries_from_spec_full_new(new_resp*spec0,factor)            
                vec_tseries[:,2] = tseries.transpose()

                resptt = resp_data3[:,kk].transpose()
                new_resp = expand_resp_new(resptt)
                tseries = synth_tseries_from_spec_full_new(new_resp*spec0,factor)            
                vec_tseries[:,3] = tseries.transpose()
                
                #vlim[int(jj/(1+doma_dims[0])),:] = np.max(np.abs(vec_tseries[:,3]),axis=0)
                #vlim[int(jj),:] = np.max(np.abs(vec_tseries[:,3]),axis=0)
                
            vlim[jj,:] = vec_tseries[:,3]
            rlim[jj,:] = vec_tseries[:,2]
            if(((jj+1)==resp_dims[0])):
                print(jj)
            
            '''
            vlim[jj,:] = vec_tseries[:,3]
            '''
            
        #+print(int(jj/(1+int(np.sqrt(resp_dims[0])))))
            #domat = doma[:,:,jj*8]
            #print(int(jj/(1+doma_dims[0])))
            #if not (jj%int(resp_dims[0]/1)):
            if 0:
                
                index = np.ravel_multi_index((jj,0),doma_dims)
                doma_data = np.array(doma_dataset[doma_list_keys[index]]['value']).view(complex)
                
                index = np.ravel_multi_index((jj,1),doma_dims)
                doma_data0 = np.array(doma_dataset[doma_list_keys[index]]['value']).view(complex)
                
                index = np.ravel_multi_index((jj,2),doma_dims)
                doma_data1 = np.array(doma_dataset[doma_list_keys[index]]['value']).view(complex)
                
                doma_data2 = doma_data0 + doma_data1 + doma_data
                
                            
                for kk in range(int(doma_data.shape[1])):
                    
                    domatt = doma_data[:,kk].transpose()
                    
                    new_doma = expand_resp(domatt)*spec0
                    tseries = synth_tseries_from_spec_full_new(new_doma,factor)
                    mat_tseries[:,kk,0] = tseries.transpose()
                    domatt = doma_data0[:,kk].transpose()
                    new_doma = expand_resp(domatt)*spec0
                    tseries = synth_tseries_from_spec_full_new(new_doma,factor)            
                    mat_tseries[:,kk,1] = tseries.transpose()
                    domatt = doma_data1[:,kk].transpose()
                    new_doma = expand_resp(domatt)*spec0
                    tseries = synth_tseries_from_spec_full_new(new_doma,factor)            
                    mat_tseries[:,kk,2] = tseries.transpose()
                    domatt = doma_data2[:,kk].transpose()
                    new_doma = expand_resp(domatt)*spec0
                    tseries = synth_tseries_from_spec_full_new(new_doma,factor)            
                    mat_tseries[:,kk,3] = tseries.transpose()             
                mlim = np.amax(np.abs(mat_tseries))/4
                for ll in range(int(factor*int(spec_size/dtsr))):
                    
                    this_series = mat_tseries[(ll+1)*dtsr-1,:,0].real/mlim
                    #this_series = 20*np.log10(abs(mat_tseries[ll*dtsr,:,0]))
                    this_series.shape = this_grid.shape
                    h0.set_array(this_series.ravel())
                    
                    this_series = mat_tseries[(ll+1)*dtsr-1,:,1].real/mlim
                    #this_series = 20*np.log10(abs(mat_tseries[ll*dtsr,:,1]))
                    this_series.shape = this_grid.shape            
                    h1.set_array(this_series.ravel())
                    
                    this_series = mat_tseries[(ll+1)*dtsr-1,:,2].real/mlim
                    #this_series = 20*np.log10(abs(mat_tseries[ll*dtsr,:,2]))
                    this_series.shape = this_grid.shape            
                    h2.set_array(this_series.ravel())
                    
                    this_series = mat_tseries[(ll+1)*dtsr-1,:,3].real/mlim
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
                    #plt.grid(True)
                    fg.canvas.flush_events()
        
    spec_grid,ndx_grid = np.meshgrid(x_spec_full,range(resp_dims[0]))
    #ax.pcolormesh(spec_grid,ndx_grid,np.abs(vlim),shading='nearest',vmin=-1, vmax=1)
    ax.pcolormesh(spec_grid,ndx_grid,rlim.real,shading='nearest')
    #ax.imshow(vlim.real)
    fg,ax = plt.subplots(1,1,figsize=(10,7))
    ax.pcolormesh(spec_grid,ndx_grid,vlim.real,shading='nearest')
    fg,ax = plt.subplots(1,1,figsize=(10,7))
    
    rng = range(3,doma_dims[0],40)
    xx, yy = np.meshgrid(rng,rng)
    rngstacked = np.vstack((xx.ravel(),yy.ravel()))
    multi_index = np.ravel_multi_index(rngstacked,(doma_dims[0],doma_dims[0]))
    
    #ax.plot(np.abs(vlim[:doma_dims[0]:]))
    vlimmax = np.max(np.abs(vlim),axis=-1)

    
    
    '''
    vlimmax.shape = (doma_dims[0],doma_dims[0])
    '''
    
    

    plt.ioff()
    #ax.pcolormesh(np.abs(vlim.transpose()))
    ax.plot(vlimmax/np.max(vlimmax))
    fg.suptitle(filename[:-3])
    x_table = np.arange(-int(resp_dims[0])/2,int(resp_dims[0])/2)[:]*0.745
    y_table = vlimmax[:]/np.max(vlimmax)
    d_table = np.array([x_table,y_table]).transpose()
    print(x_table.shape)
    print(y_table.shape)
    np.savetxt('vlimmax.csv',d_table, delimiter=',', header=','.join(('t','s')), comments='')
    plt.show()
    fg,ax = plt.subplots(1,1,figsize=(10,7))
    
    #plt.ioff()
    spec_grid,ndx_grid = np.meshgrid(x_spec_full,range(256))
    #ax.pcolormesh(spec_grid,ndx_grid,np.max(np.abs(vlim,axis=0),shading='nearest',vmin=-1, vmax=1)
    '''
    ax.pcolormesh(spec_grid,ndx_grid,np.abs(vlim[multi_index,:]),shading='nearest',vmin=-1, vmax=1)
    plt.show()
    '''