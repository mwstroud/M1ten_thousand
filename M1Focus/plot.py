from bmtk.analyzer.compartment import plot_traces
from bmtk.analyzer.spike_trains import plot_raster
from bmtk.analyzer.spike_trains import spike_statistics
from bmtk.analyzer.spike_trains import plot_rates_boxplot
from scipy.signal import hanning,welch,decimate
import h5py
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import sys
import math
from bmtools.cli.plugins.util import util
from os.path import exists
from bmtk.analyzer.spike_trains import plot_raster
import os, sys

def raster(spikes_df,node_set,skip_ms=0,ax=None):
    spikes_df = spikes_df[spikes_df['timestamps']>skip_ms] 
    for node in node_set:
        cells = node['ids'] 
        
        cell_spikes = spikes_df[spikes_df['node_ids'].isin(cells)]
        #print(cell_spikes)
        ax.scatter(cell_spikes['timestamps'],cell_spikes['node_ids'],
                   c='tab:'+node['color'],s=0.25, label=node['name'])
    
    handles,labels = ax.get_legend_handles_labels()
    ax.set_title('M1 Cortex Raster')
    ax.legend(reversed(handles), reversed(labels))
    ax.grid(True)

def spike_frequency_histogram(spikes_df,node_set,ms,skip_ms=0,ax=None,n_bins=10):
    print("Type : mean (std)")
    for node in node_set:
        cells = node['ids'] 
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
            ax.set_title('M1 Cortex Rates')
    if ax:
        ax.set_xscale('log')
        ax.legend() 

def populations(config):
    nodes = util.load_nodes_from_config(config)
    i=0
    j=0
    CP_nodes = []
    CS_nodes = []
    FSI_nodes = []
    LTS_nodes = []
    for j in nodes:
        node=nodes[j]
        num_cells=len(node['node_type_id'])
        if node['model_type'][0]=='biophysical':
            #print(node)
            #print(node['node_id'])
            num_cells=len(node['node_type_id'])
            for i in range(num_cells-1):
                if(node['pop_name'][i]=='CP'):
                    CP_nodes.append(i)
                if(node['pop_name'][i]=='CS'):
                    CS_nodes.append(i)
                if(node['pop_name'][i]=='FSI'):
                    FSI_nodes.append(i)
                if(node['pop_name'][i]=='LTS'):
                    LTS_nodes.append(i)
    return CP_nodes, CS_nodes, FSI_nodes, LTS_nodes 

def core_populations(config):
    x_core=[250,350]
    y_core=[250,350]
    z_core=[300,1500]
    nodes = util.load_nodes_from_config(config)
    i=0
    j=0
    CP_nodes = []
    CS_nodes = []
    FSI_nodes = []
    LTS_nodes = []
    for j in nodes:
        node=nodes[j]
        num_cells=len(node['node_type_id'])
        if node['model_type'][0]=='biophysical':
            #print(node)
            #print(node['node_id'])
            num_cells=len(node['node_type_id'])
            core = False
            for i in range(num_cells-1):
                if node['pos_x'][i] > x_core[0] and node['pos_x'][i] < x_core[1] and node['pos_y'][i] > y_core[0] and node['pos_y'][i] < y_core[1] and node['pos_z'][i] > z_core[0] and node['pos_z'][i] < z_core[1]:
                    core = True
                if(node['pop_name'][i]=='CP'):
                    if core == True:
                        CP_nodes.append(i)
                if(node['pop_name'][i]=='CS'):
                    if core == True:
                        CS_nodes.append(i)
                if(node['pop_name'][i]=='FSI'):
                    if core == True:
                        FSI_nodes.append(i)
                if(node['pop_name'][i]=='LTS'):
                    if core == True:
                        LTS_nodes.append(i)

    return CP_nodes, CS_nodes, FSI_nodes, LTS_nodes 

def intermediate_populations(config):
    x_core=[250,350]
    y_core=[250,350]
    z_core=[300,1500]
    nodes = util.load_nodes_from_config(config)
    i=0
    j=0
    CP_nodes = []
    CS_nodes = []
    FSI_nodes = []
    LTS_nodes = []
    for j in nodes:
        node=nodes[j]
        num_cells=len(node['node_type_id'])
        if node['model_type'][0]=='biophysical':
            #print(node)
            #print(node['node_id'])
            num_cells=len(node['node_type_id'])
            core = False
            for i in range(num_cells-1):
                if node['pos_x'][i] > x_core[0] and node['pos_x'][i] < x_core[1] and node['pos_y'][i] > y_core[0] and node['pos_y'][i] < y_core[1] and node['pos_z'][i] > z_core[0] and node['pos_z'][i] < z_core[1]:
                    core = False
                if(node['pop_name'][i]=='CP'):
                    if core == False:
                        CP_nodes.append(i)
                if(node['pop_name'][i]=='CS'):
                    if core == False:
                        CS_nodes.append(i)
                if(node['pop_name'][i]=='FSI'):
                    if core == False:
                        FSI_nodes.append(i)
                if(node['pop_name'][i]=='LTS'):
                    if core == False:
                        LTS_nodes.append(i)
    return CP_nodes, CS_nodes, FSI_nodes, LTS_nodes
                    
def gaussian(x,mean,std,pmax):
    scale= pmax*std*math.sqrt(2*math.pi)
    y=(scale/(std*math.sqrt(2*math.pi)))*math.exp(-0.5*((x-mean)**2/(2*std**2)))   
    return y 
                

def plot(choose):
    #CS_nodes=list(range(0,10))+list(range(128,148))+list(range(276,386))+list(range(619,709))
    #CP_nodes=list(range(266,286))+list(range(429,628))
    #FSI_nodes=list(range(88,122))+list(range(226,260))+list(range(389,422))+list(range(720,754))+list(range(950,984))
    #LTS_nodes=list(range(112,138))+list(range(250,276))+list(range(413,439))+list(range(744,769))+list(range(974,1000))
    
    #CC_nodes=list(range(10,98))+list(range(138,236))+list(range(381,398))+list(range(704,729))
    #CTH_nodes=list(range(376,391))+list(range(699,714))+list(range(760,960))



    if choose=='1':
        print("plotting spike raster")
        [CP_nodes, CS_nodes, FSI_nodes, LTS_nodes] = populations("config.json")
        node_set = [
            {"name":"CP","ids":CP_nodes,"color":"blue"},
            {"name":"CS","ids":CS_nodes,"color":"red"},
            {"name":"FSI","ids":FSI_nodes,"color":"green"},
            {"name":"LTS","ids":LTS_nodes,"color":"purple"}
        ]

        if exists('output/spikes_short.h5'):
            spikes_location = 'output/spikes_short.h5'
        elif exists('output/spikes_long.h5'):
            spikes_location = 'output/spikes_long.h5'
        elif exists('output/spikes_baseline.h5'):
            spikes_location = 'output/spikes_baseline.h5'
        f = h5py.File(spikes_location)
        spikes_df = pd.DataFrame({'node_ids':f['spikes']['cortex']['node_ids'],'timestamps':f['spikes']['cortex']['timestamps']}) 

        skip_ms=0
        end_ms = 4000
        fig, ax = plt.subplots(1,1,figsize=(6,5))#6.4,4.8 default
        #_ = plot_raster(config_file='config.json',group_by='pop_name',group_excludes=['CC','CTH'])
        raster(spikes_df,node_set,ax=ax, skip_ms=skip_ms)
        plt.show()
    elif choose=='2':
        print("plotting firing rates")
        #[CP_nodes, CS_nodes, FSI_nodes, LTS_nodes] = populations("config.json")
        #node_set_total = [
            #{"name":"CP","ids":CP_nodes,"color":"blue"},
            #{"name":"CS","ids":CS_nodes,"color":"red"},
            #{"name":"FSI","ids":FSI_nodes,"color":"green"},
            #{"name":"LTS","ids":LTS_nodes,"color":"purple"}
        #]
        [cCP_nodes, cCS_nodes, cFSI_nodes, cLTS_nodes] = core_populations("config.json")
        node_set_core = [
            {"name":"CP","ids":cCP_nodes,"color":"blue"},
            {"name":"CS","ids":cCS_nodes,"color":"red"},
            {"name":"FSI","ids":cFSI_nodes,"color":"green"},
            {"name":"LTS","ids":cLTS_nodes,"color":"purple"}
        ]
        [iCP_nodes, iCS_nodes, iFSI_nodes, iLTS_nodes] = intermediate_populations("config.json")
        node_set_intermediate = [
            {"name":"CP","ids":iCP_nodes,"color":"blue"},
            {"name":"CS","ids":iCS_nodes,"color":"red"},
            {"name":"FSI","ids":iFSI_nodes,"color":"green"},
            {"name":"LTS","ids":iLTS_nodes,"color":"purple"}
        ]
        

        if exists('output/spikes_short.h5'):
            spikes_location = 'output/spikes_short.h5'
        elif exists('output/spikes_long.h5'):
            spikes_location = 'output/spikes_long.h5'
        elif exists('output/spikes_baseline.h5'):
            spikes_location = 'output/spikes_baseline.h5'
        f = h5py.File(spikes_location)
        spikes_df = pd.DataFrame({'node_ids':f['spikes']['cortex']['node_ids'],'timestamps':f['spikes']['cortex']['timestamps']}) 

        skip_ms=0
        end_ms = 4000
        fig, ax = plt.subplots(1,2,figsize=(12,5))#6.4,4.8 default
        #_ = plot_rates_boxplot(config_file='config.json',group_by='pop_name',group_excludes=['CC','CTH'])
        #spike_frequency_histogram(spikes_df,node_set_total,end_ms,ax=ax[0],skip_ms=skip_ms)
        spike_frequency_histogram(spikes_df,node_set_core,end_ms,ax=ax[0],skip_ms=skip_ms)
        ax[0].set_title("Core Rates")
        spike_frequency_histogram(spikes_df,node_set_intermediate,end_ms,ax=ax[1],skip_ms=skip_ms)
        ax[1].set_title("Intermediate Rates")
        plt.show()
    elif choose=='3':
        print("plotting voltage trace")
        _ = plot_traces(config_file='config.json', report_name='v_report',group_by='pop_name',group_excludes=['CP','CS','LTS','CC','CTH'], average=True)
    elif choose=='4':
        print("plotting syn report")
        _ = plot_traces(report_path='./output/Gfluct_exc.h5')
        #./output/PN2FSI_nmda_report.h5
        #./output/PN2FSI_ampa_report.h5
        #./output/Gfluct_exc.h5
        #./output/Gfluct_inh.h5
    elif choose=='5':
        x=range(-300, 300)
        gaus= np.vectorize(gaussian)
        y=gaus(x,0,124.62,0.05)
        plt.plot(x,y)
        plt.show()
    else:
        return

if __name__ == '__main__':
    if sys.argv[1]:
        plot(choose=sys.argv[1])
        #intermediate_populations('config.json')
    else:
        #intermediate_populations('config.json')
        plot(choose=1)