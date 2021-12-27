import numpy as np
import sys
from bmtools.cli.plugins.util import util

def populations(config):
    nodes = util.load_nodes_from_config(config)
    CP_nodes = []
    CS_nodes = []
    FSI_nodes = []
    LTS_nodes = []
    for j in nodes:
        node=nodes[j]
        num_cells=len(node['node_type_id'])
        if node['model_type'][0]=='virtual':
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

def print_cell_ranges(config):
    CP_nodes, CS_nodes, FSI_nodes, LTS_nodes = populations(config);
    print("CP_nodes:")
    print(CP_nodes)
    print("CS_nodes:")
    print(CS_nodes)
    print("FSI_nodes:")
    print(FSI_nodes)
    print("LTS_nodes:")
    print(LTS_nodes)

print_cell_ranges('config.json')
