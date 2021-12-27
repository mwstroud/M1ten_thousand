import sys
from bmtools.cli.plugins.util import util
from bmtk.utils.reports.spike_trains import PoissonSpikeGenerator
import random
import csv
import numpy as np
import time
from console_progressbar import ProgressBar


#################################################################
####### IMPORTANT: the code here doesnt scale properly! #########
#################################################################

def lognorm_fr_list(n,m,s):
    mean = np.log(m) - 0.5 * np.log((s/m)**2+1)
    std = np.sqrt(np.log((s/m)**2 + 1))
    ranlog= np.random.lognormal(mean,std)
    #print(ranlog)
    return [ranlog for i in range(n)]

#def build_input(t_sim, numPN_A = 640, numPN_C=260, numBask = 100):
def build_poisson_input(population,node_ids,mean,std,output_h5,tstart,tend,t_sim=15000):
    print('Building input for ' + population + "[" + str(len(node_ids)) + " cells at " + str(mean) + "(" + str(std) + ") Hz]")
    psg = PoissonSpikeGenerator(population=population)
    psg.add(node_ids=node_ids,  
    firing_rate=lognorm_fr_list(len(node_ids),mean,std),
    times=(tstart, tend))  
    psg.to_sonata(output_h5)

def populations(config):
    nodes = util.load_nodes_from_config(config)
    CP_nodes = []
    CS_nodes = []
    FSI_nodes = []
    LTS_nodes = []
    for j in nodes:
        node=nodes[j]
        num_cells=len(node['node_type_id'])
        if node['model_type'][0]=='biophysical':
            num_cells=len(node['node_type_id'])
            CP_id = []
            CS_id = []
            FSI_id = []
            LTS_id = []
            for i in range(num_cells):
                if(node['pop_name'][i]=='CP'):
                    CP_id.append(i)
                    CP_nodes.append(node.iloc[i])
                if(node['pop_name'][i]=='CS'):
                    CS_id.append(i)
                    CS_nodes.append(node.iloc[i])
                if(node['pop_name'][i]=='FSI'):
                    FSI_id.append(i)
                    FSI_nodes.append(node.iloc[i])
                if(node['pop_name'][i]=='LTS'):
                    LTS_id.append(i)
                    LTS_nodes.append(node.iloc[i])
    #print(len(CP_nodes))
    #print(len(CS_nodes))
    #print(len(FSI_nodes))
    #print(len(LTS_nodes))
    
    return CP_nodes, CS_nodes, FSI_nodes, LTS_nodes, CP_id, CS_id, FSI_id, LTS_id

def thalnodes(config):
    nodes = util.load_nodes_from_config(config)
    Thal_nodes = []
    for j in nodes:
        node=nodes[j]
        num_cells=len(node['node_type_id'])
        if node['model_type'][0]=='virtual':
            num_cells=len(node['node_type_id'])
            for i in range(num_cells):
                if(node['pop_name'][i]=='thal'):
                    Thal_nodes.append(i)
    return Thal_nodes 

def Intnodes(config):
    nodes = util.load_nodes_from_config(config)
    Int_nodes = []
    for j in nodes:
        node=nodes[j]
        num_cells=len(node['node_type_id'])
        if node['model_type'][0]=='virtual':
            num_cells=len(node['node_type_id'])
            for i in range(num_cells):
                if(node['pop_name'][i]=='Int'):
                    Int_nodes.append(i)
    return Int_nodes


def coreChoice(virt,CP_nodes,CS_nodes,CP_id,CS_id):
    if len(virt) == len(CP_nodes)+len(CS_nodes):
        x1 = np.array([CP_nodes[i]['pos_x'] for i in range(len(CP_nodes))])
        y1 = np.array([CP_nodes[i]['pos_y'] for i in range(len(CP_nodes))])
        z1 = np.array([CP_nodes[i]['pos_z'] for i in range(len(CP_nodes))])
        
        x2 = np.array([CS_nodes[i]['pos_x'] for i in range(len(CS_nodes))])
        y2 = np.array([CS_nodes[i]['pos_y'] for i in range(len(CS_nodes))])
        z2 = np.array([CS_nodes[i]['pos_z'] for i in range(len(CS_nodes))])
        coreCP = []
        coreCS = []
        # These for loops will go through the CS and CP nodes and generate a list of corresponding virtual cell ids that have PNs in the core
        tempCP = []
        tempCS = []
        for i in range(len(CP_nodes)):
            if x1[i] >=300 and x1[i] <=700 and y1[i] >=300 and y1[i] <=700 and z1[i] >=300 and z1[i] <=700:
                #print('CP node position: ' + str(CP_nodes[i]['pos_x']) + ' , ' + str(CP_nodes[i]['pos_y']) + ' , ' + str(CP_nodes[i]['pos_z']))
                #print(virt[i])
                if CP_id[i] > 4000:
                    CP_id[i] -= 1000
                tempCP.append(CP_id[i])
        for i in range(len(CS_nodes)):
            if x2[i] >=300 and x2[i] <=700 and y2[i] >=300 and y2[i] <=700 and z2[i] >=300 and z2[i] <=700:
                #print('CS node position: ' + str(CS_nodes[i]['pos_x']) + ' , ' + str(CS_nodes[i]['pos_y']) + ' , ' + str(CS_nodes[i]['pos_z']))
                if CS_id[i] > 4000:
                    CS_id[i] -= 1000
                tempCS.append(CS_id[i])
        #print(tempCP)
        # I got the cell ids of the CP and CS cells that are in the core, now to match it to the virtual cell ids
        for i in range(len(virt)):
            if np.isin(i,tempCP):
                coreCP.append(i)
            elif np.isin(i,tempCS):
                coreCS.append(i)
                
        # Now I just need to divide the coreids by 8 and return a
        assemblies = 8
        num_coreCP = len(coreCP)
        num_coreCS = len(coreCS)
        #print('# of Core CP cells: '+ str(num_coreCP))
        #print('# of Core CS cells: '+ str(num_coreCS))
        
        # Shuffle the nodes
        #random.shuffle(Thal_nodes)
        random.shuffle(coreCP)
        random.shuffle(coreCS)
        #print(Thal_nodes)
        thal1=[]
        thal2=[]
        thal3=[]
        thal4=[]
        thal5=[]
        thal6=[]
        thal7=[]
        thal8=[]
        # Split the CP nodes into 8 assemblies
        for i in range(num_coreCP):
            if i < round(1*(num_coreCP/assemblies)):
                thal1.append(coreCP[i])
            if i >= round(1*(num_coreCP/assemblies)) and i < round(2*(num_coreCP/assemblies)):
                thal2.append(coreCP[i])
            if i >= round(2*(num_coreCP/assemblies)) and i < round(3*(num_coreCP/assemblies)):
                thal3.append(coreCP[i])
            if i >= round(3*(num_coreCP/assemblies)) and i < round(4*(num_coreCP/assemblies)):
                thal4.append(coreCP[i])
            if i >= round(4*(num_coreCP/assemblies)) and i < round(5*(num_coreCP/assemblies)):
                thal5.append(coreCP[i])
            if i >= round(5*(num_coreCP/assemblies)) and i < round(6*(num_coreCP/assemblies)):
                thal6.append(coreCP[i])
            if i >= round(6*(num_coreCP/assemblies)) and i < round(7*(num_coreCP/assemblies)):
                thal7.append(coreCP[i])
            if i >= round(7*(num_coreCP/assemblies)) and i < round(8*(num_coreCP/assemblies)):
                thal8.append(coreCP[i])
                
        # Split the CS nodes into 8 assemblies. Append to the same thal lists 
        for i in range(num_coreCS):
            if i < round(1*(num_coreCS/assemblies)):
                thal1.append(coreCS[i])
            if i >= round(1*(num_coreCS/assemblies)) and i < round(2*(num_coreCS/assemblies)):
                thal2.append(coreCS[i])
            if i >= round(2*(num_coreCS/assemblies)) and i < round(3*(num_coreCS/assemblies)):
                thal3.append(coreCS[i])
            if i >= round(3*(num_coreCS/assemblies)) and i < round(4*(num_coreCS/assemblies)):
                thal4.append(coreCS[i])
            if i >= round(4*(num_coreCS/assemblies)) and i < round(5*(num_coreCS/assemblies)):
                thal5.append(coreCS[i])
            if i >= round(5*(num_coreCS/assemblies)) and i < round(6*(num_coreCS/assemblies)):
                thal6.append(coreCS[i])
            if i >= round(6*(num_coreCS/assemblies)) and i < round(7*(num_coreCS/assemblies)):
                thal7.append(coreCS[i])
            if i >= round(7*(num_coreCS/assemblies)) and i < round(8*(num_coreCS/assemblies)):
                thal8.append(coreCS[i])


        
        # Short burst input for 1000 ms with 100 ms subbursts followed by 500 ms silence
        Thal=[thal1,thal2,thal3,thal4,thal5,thal6,thal7,thal8]
        return Thal
    else:
        print("Number of Thalamus cells don't match number of PNs")
        return -1
     

def build_input(t_sim, num_thal = 800, num_CS=200, num_CTH = 200, num_CC=200, num_FSI=120, num_LTS=80, scale=1):
    pb = ProgressBar(total=100, decimals=3, length=50, fill='X', zfill='-')
    print("Building all input spike trains")
    pb.print_progress_bar(2)
    
    # Get nodes                   
    Int_nodes=Intnodes("config.json")
    num_int=len(Int_nodes)
    Thal_nodes=thalnodes("config.json")
    num_thal = len(Thal_nodes)
    CP_nodes, CS_nodes, FSI_nodes, LTS_nodes, CP_id, CS_id, FSI_id, LTS_id = populations("config.json")
    #print(CP_id)
    pb.print_progress_bar(10)
    #print(CP_nodes)
    # Get the thalamic node id that will connect to CP vs CS cells
    # CP_nodes = []
    # CS_nodes = []
    # for i in range(num_thal):
    #     if  Thal_nodes[i] < 200:
    #         CP_nodes.append(Thal_nodes[i])
    #     elif Thal_nodes[i] >= 200 and Thal_nodes[i] < 4000:
    #         CS_nodes.append(Thal_nodes[i])
    #     elif Thal_nodes[i] >= 4000 and Thal_nodes[i] < 7800:
    #         CP_nodes.append(Thal_nodes[i])
    #     elif Thal_nodes[i] >= 7800 and Thal_nodes[i] < 8000:
    #         CS_nodes.append(Thal_nodes[i])
    
    num_CP = len(CP_nodes)
    num_CS = len(CS_nodes)
    assemblies = 8
    num_group = round(num_thal/assemblies)
    #print(len(Thal_nodes))
    #print(len(CP_nodes))
    Thal = coreChoice(Thal_nodes, CP_nodes, CS_nodes, CP_id, CS_id)
    num_chosen = len(Thal[0]) + len(Thal[1]) + len(Thal[2]) + len(Thal[3]) + len(Thal[4]) + len(Thal[5]) + len(Thal[6]) + len(Thal[7])
    pb.print_progress_bar(20)
    assemblies = 8
    with open("./input/Assembly_ids.csv", "w", newline="") as f:
        writer = csv.writer(f)
        for row in Thal:
            writer.writerow(row)
        
    psgs = PoissonSpikeGenerator(population='thalamus')
    for i in range(8):
        for j in range(assemblies):
            time1=i*1.5 + j*0.125
            time2=i*1.5 + j*0.125 +0.125
            psgs.add(node_ids=Thal[j],  
            firing_rate=lognorm_fr_list(num_group*scale,50,0),
            times=(time1, time2)) 
            pb.print_progress_bar(30+2*i)
    # Long burst input for 1000 ms followed by 500 ms silence
    
    psgl = PoissonSpikeGenerator(population='thalamus') 
    for i in range(assemblies):
        time1=i*1.5
        time2=i*1.5+1
        psgl.add(node_ids=Thal[i],  
            firing_rate=lognorm_fr_list(num_group*scale,50,0),
            times=(time1, time2))
        pb.print_progress_bar(50+2*i)   
    #print(len(Thal_nodes))
    # Thalamus test constant 50 Hz input
    psgc = PoissonSpikeGenerator(population='thalamus')
    psgc.add(node_ids=Thal_nodes,  
    firing_rate=lognorm_fr_list(num_thal*scale,50,0),
    times=(0, 10))
    pb.print_progress_bar(70)  
    
    # Thalamus baseline
    psgb = PoissonSpikeGenerator(population='thalamus')
    psgb.add(node_ids=Thal_nodes,  
    firing_rate=lognorm_fr_list(num_thal*scale,2,0),
    times=(0, 10))  
    pb.print_progress_bar(80)
    # Interneuron baseline
    psgi = PoissonSpikeGenerator(population='Intbase')
    psgi.add(node_ids=Int_nodes,  
    firing_rate=lognorm_fr_list(num_int*scale,2,0),
    times=(0, 10))
    pb.print_progress_bar(90)
    
    psgc.to_sonata('./input/thalamus_const.h5')
    psgb.to_sonata('./input/thalamus_base.h5')
    psgs.to_sonata('./input/thalamus_short.h5')
    psgl.to_sonata('./input/thalamus_long.h5')
    psgi.to_sonata('./input/Intbase.h5')
    # These inputs are for the baseline firing rates of the cells in the shell.
    ############ Short Burst Input ##############
    num_CP = len(CP_id)
    num_CS = len(CS_id)
    num_FSI = len(FSI_id)
    num_LTS = len(LTS_id)
    # CP Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=CP_id,  
    firing_rate=lognorm_fr_list(num_CP*scale,0.5,0.2),
    times=(0, 11.5))  
    psgl.to_sonata('./input/CP_shell_short.h5')

    # CS Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=CS_id,  
    firing_rate=lognorm_fr_list(num_CS*scale,0.5,0.2),
    times=(0, 11.5))  
    psgl.to_sonata('./input/CS_shell_short.h5')

    # FSI Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=FSI_id,  
    firing_rate=lognorm_fr_list(num_FSI*scale,5,3),
    times=(0, 11.5))
    psgl.to_sonata('./input/FSI_shell_short.h5')

    # LTS Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=LTS_id,  
    firing_rate=lognorm_fr_list(num_LTS*scale,1,0.1),
    times=(0, 11.5))
    psgl.to_sonata('./input/LTS_shell_short.h5')
    ######### Long Burst Shell input ##########
    # CP Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=CP_id,  
    firing_rate=lognorm_fr_list(num_CP*scale,0.5,0.2),
    times=(0, 11.5))  
    psgl.to_sonata('./input/CP_shell_long.h5')

    # CS Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=CS_id,  
    firing_rate=lognorm_fr_list(num_CS*scale,0.5,0.2),
    times=(0, 11.5))  
    psgl.to_sonata('./input/CS_shell_long.h5')

    # FSI Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=FSI_id,  
    firing_rate=lognorm_fr_list(num_FSI*scale,5,3),
    times=(0, 11.5))
    psgl.to_sonata('./input/FSI_shell_long.h5')

    # LTS Average
    psgl = PoissonSpikeGenerator(population='shell')
    psgl.add(node_ids=LTS_id,  
    firing_rate=lognorm_fr_list(num_LTS*scale,1.5,1),
    times=(0, 11.5))
    psgl.to_sonata('./input/LTS_shell_long.h5')
    pb.print_progress_bar(100)
    print('Number of Thal nodes activated: '+ str(num_chosen))
    print("Done")


if __name__ == '__main__':
    if __file__ != sys.argv[-1]:
        build_input(int(sys.argv[-1]))
    else:
        build_input(10)
