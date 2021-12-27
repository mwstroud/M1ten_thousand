"""Script for running the network built in build_network.py
Also saves a file called Connections.csv that consists of information about
each synapse in the simulation.
"""

import os, sys, inspect

currentdir = os.path.dirname(os.path.abspath(inspect.getfile(inspect.currentframe())))
parentdir = os.path.dirname(currentdir)
sys.path.insert(0, parentdir)
sys.path.insert(0, currentdir)

from bmtk.simulator import bionet
import numpy as np
from neuron import h
import synapses
import pandas as pd



def save_connections(graph, sim):
    """Saves Connections.csv based on the given network.
    Parameters
    ----------
    graph : BioNetwork
        the network that the connections are retrieved from
    sim : BioSimulator
        the simulation about to be run (not used in this function)
    """
    cells = graph.get_local_cells()
    cell = cells[list(cells.keys())[0]]

    h.distance(sec=cell.hobj.soma[0])  # Makes the distances correct.

    sec_types = []  # soma, apic, or dend
    weights = []  # scaled conductances (initW)
    dists = []  # distance from soma
    node_ids = []  # node id within respective node population (exc, prox_inh, dist_inh)
    names = []  # full NEURON str representation of postsynaptic segment
    source_pops = []  # node population
    release_probs = []  # propability of release.

    for c in cell.connections():
        con = c._connector
        source = c.source_node
        syn = con.syn()
        seg = con.postseg()
        print (seg.keys())
        fullsecname = seg.sec.name()
        #source_pops.append(source.node_gid)
        #node_ids.append(source._node_id)
        #weights.append(float(syn.initW))
        #release_probs.append(float(syn.P_0))
        names.append(str(seg))
        sec_types.append(fullsecname.split(".")[1][:4])
        dists.append(float(h.distance(seg)))

    df = pd.DataFrame()
    #df["Node ID"] = node_ids
    #df["Conductance"] = weights
    df["Type"] = sec_types
    df["Name"] = names
    #df["SecID"] = 
    #df["Source Population"] = source_pops
    #df["Release Probability"] = release_probs
    df.to_csv("Connections.csv", index=False)

synapses.load()
config_file = 'config.json'
conf = bionet.Config.from_json(config_file, validate=True)
conf.build_env()
graph = bionet.BioNetwork.from_config(conf)
sim = bionet.BioSimulator.from_config(conf, network=graph)

save_connections(graph,sim)