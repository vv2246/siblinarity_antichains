# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 12:52:55 2018
@author: Vaiva & Tim
"""


from collections import defaultdict
import itertools
import networkx as nx
import numpy as np
import scipy.special as ss
from scipy.stats import hypergeom
import math    
import random
import matplotlib.pyplot as plt
import os
import copy



def tr(DAG, output=False):
    """
    Transitive reduction,courtesy of JC
    """
    # for printing progress
    E = DAG.number_of_edges()
    i = 0
    print_limit = 10
    print_counter = print_limit
    edges = list(DAG.edges())
    #########################
    for edge in edges:
        # check edge is necessary for causal structure
        [a, b] = edge
        DAG.remove_edge(a, b)
        if not nx.has_path(DAG, a, b):
            DAG.add_edge(a, b)
        
        if output:
            i += 1
            pc = (i/float(E))*100
            if pc > print_counter:
                print ('Finished %s percent' % int(math.floor(pc)))
                print_counter += print_limit
                
    return DAG


    


def is_antichain(graph,nodelist):
    """
    Tests whether a list of nodes in a graph is an antichain.
    (TSE tidy version)
    Note this works for any type of graph but only make sense for directed or directed acyclic graphs.
    For a simple graph, the maximal antichains are just the components.
    
    Parameters
    ----------
    graph = networkx graph 
    nodelist = list of nodes suggested as an antichain
    
    Returns
    -------
    Bool: 
        True if the nodelist is an antichain in graph
        False if the nodelist is not an antichain in graph
    """
    for n,m in itertools.combinations(nodelist,2):
        if nx.has_path(graph,m,n) or nx.has_path(graph,n,m):
            return False
    return True  


def is_weakly_connected(graph,source_nodes,target_nodes):
    """
    Tests whether a list of source nodes in a graph have a path in either direction between a list of target nodes.
        
    Parameters
    ----------
    graph = networkx graph 
    source_nodes = list of source nodes for paths
    target_nodes = list of target nodes for paths
    
    Returns
    -------
    Bool: 
        True if there is a path from at least one source node to at least one target node or the other way round
        False otherwise
    """
    for s,t in itertools.product(source_nodes,target_nodes):
        if nx.has_path(graph,s,t) or nx.has_path(graph,t,s):
            return True
    return False  