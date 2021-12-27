"""
A python implementation of matlab/analysis.m

TB - 8/4/21
"""

from scipy.signal import hanning,welch,decimate
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys

from bmtk.analyzer.spike_trains import plot_raster






def raster(spikes_df,node_set,skip_ms=0,ax=None):
    spikes_df = spikes_df[spikes_df['timestamps']>skip_ms] 
    for node in node_set:
        cells = range(node['start'],node['end']+1) #+1 to be inclusive of last cell
        cell_spikes = spikes_df[spikes_df['node_ids'].isin(cells)]

        ax.scatter(cell_spikes['timestamps'],cell_spikes['node_ids'],
                   c='tab:'+node['color'],s=0.25, label=node['name'])
    
    handles,labels = ax.get_legend_handles_labels()
    ax.legend(reversed(handles), reversed(labels))
    ax.grid(True)

def raw_ecp(lfp):
    pass

def ecp_psd(ecp,skip_n=0,downsample=20,nfft=1024,fs=1000,noverlap=0,ax=None):
    
    #skip_n first few
    ecp = [x for x in ecp if np.isnan(x) == False]
    #print(ecp)
    data = ecp[skip_n:]
    #print(data)
    #downsample the data to fit ms (steps used 20=1/.05 step)
    lfp_d = decimate(ecp,downsample)
    raw_ecp(lfp_d)
    win = hanning(nfft, True)

    f,pxx = welch(lfp_d,fs,window=win,noverlap=noverlap,nfft=nfft)

    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.plot(f, pxx*1000,linewidth=0.6)
    ax.set_xlabel("frequency [Hz]")
    ax.set_ylabel("PSD [V**2/Hz]")
    ax.set_ylim([0.00001,0.1])
    
    
    

def spike_frequency_histogram(spikes_df,node_set,ms,skip_ms=0,ax=None,n_bins=10):
    print("Type : mean (std)")
    for node in node_set:
        cells = range(node['start'],node['end']+1) #+1 to be inclusive of last cell
        cell_spikes = spikes_df[spikes_df['node_ids'].isin(cells)]

        #skip the first few ms
        cell_spikes = cell_spikes[cell_spikes['timestamps']>skip_ms]
        spike_counts = cell_spikes.node_ids.value_counts()
        total_seconds = (ms-skip_ms)/1000
        spike_counts_per_second = spike_counts / total_seconds

        spikes_mean = spike_counts_per_second.mean()
        spikes_std = spike_counts_per_second.std()
        
        label = "{} : {:.2f} ({:.2f})".format(node['name'],spikes_mean,spikes_std)
        print(label)
        c = "tab:" + node['color']
        if ax:
            ax.hist(spike_counts_per_second,n_bins,density=True,histtype='bar',label=label,color=c)
    if ax:
        ax.set_xscale('log')
        ax.legend() 
        
        

def run(show_plots=False):
    

    dt = 0.05
    steps_per_ms = 1/dt
    skip_seconds = 5
    skip_ms = skip_seconds*1000
    skip_n = int(skip_ms * steps_per_ms)
    end_ms = 15000

    spikes_location = 'output/spikes.h5'
    
    #print("loading " + spikes_location)
    f = h5py.File(spikes_location)
    #print(f['spikes'])
    spikes_df = pd.DataFrame({'node_ids':f['spikes']['cortex']['node_ids'],'timestamps':f['spikes']['cortex']['timestamps']}) 
    #print(spikes_df)
    #print("done")

    if show_plots:
        ecp_h5_location = 'output/ecp.h5'
        print("loading " + ecp_h5_location)
        ecp_channel = 0
        f = h5py.File(ecp_h5_location)
        data_raw = np.array(f['ecp']['data'])
        ecp = data_raw.T[ecp_channel] #flip verts and grab channel 0
        print("done")

    node_set = [
        {"name":"CP","start":0,"end":799,"color":"blue"},
        {"name":"CS","start":800,"end":892,"color":"red"},
        {"name":"FSI","start":893,"end":943,"color":"green"},
        {"name":"LTS","start":944,"end":999,"color":"purple"}
    ]
    
    if show_plots:
        print("plotting...")
        fig, ax2 = plt.subplots(1,1,figsize=(6,5))#6.4,4.8 default
        fig.suptitle('M1 Cortex Analysis')
        ecp_psd(ecp, skip_n=skip_n, ax=ax2)
        #ecp = [x for x in ecp if np.isnan(x) == False]
        #plt.psd(ecp)
        spike_frequency_histogram(spikes_df,node_set,end_ms,skip_ms=skip_ms,ax=ax3)
        raster(spikes_df,node_set,skip_ms=skip_ms,ax=ax1)
        #plot_raster(config_file='simulation_config.json',group_by='pop_name')
        print("showing plots...")
        fig.tight_layout()
        plt.show()
    else:
        #spike_frequency_histogram(spikes_df,node_set,end_ms,skip_ms=skip_ms)
        pass


if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        run(False)
    else:
        run(True)
