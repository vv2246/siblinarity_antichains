# -*- coding: utf-8 -*-
"""
Created on Wed Jun 20 12:52:55 2018
Used to define models with various static methods.
@author: Vaiva
"""

#from collections import defaultdict
#import itertools
import networkx as nx
#import numpy as np
#import scipy.special as ss
import math
import random
from utilities_antichain import tr

class Model:
    # Here go class variables shared by all instances
    
    def __init__(self):
        #define instance variables here using self.<variable name>=None or whatever
        pass
    
    @staticmethod
    def f(y_1,y_2,x_1,x_2,c,N):
        """
        Calculates Manhattan distance. 
        
        Input
        y_1 : y coordinate of first point
        y_2 : y coordinate of second point
        x_1 : x coordinate of first point
        x_2 : x coordinate of second point
        c,N
        
        """
        if c<math.inf:
            norm = c/(N*(1+c))
            val = abs(y_1-y_2)+ (abs(x_1-x_2)/c)
        else:
            norm = 1/N
            val = abs(y_1-y_2)
        return norm*val
    
    @staticmethod
    def spatial_antichain_dag(N,c=1,TR=True,time_label='t',space_label='x'):    
        """
        Model of a DAG, consisting of N nodes. 
        Nodes are positioned uniformally in index order (scaled to be 0 to 1) in a time direction but
        their x coordinate is set at random between 0 and 1.
        The edges are set with probability c but only if 
        Edge probability is proportional to Manhattan distance between two nodes.
        
        Input
        N : Number of Nodes
        c : if = inf - there is no x dependence in the probability within the random model and we have only y-dependence. 
            If = 1 - the dependence on x and y are equal
            If < 1 - the dependence on y is smaller than that of x
            If > 1 - reverse
        TR=True : selects transitively reduced version if True
        time_label='t': label nodes with time coordinate
        space_label='x': label nodes with space coordinate
        
        Return
        networkx directed graph representing Tim's example DAG
        
        """
        pos= {}
        G = nx.DiGraph()
        for t in range(N):
            G.add_node(t)
            x,y = random.random(),t/N
            pos[t] = (x,y)
            G.node[t][space_label] = x
            G.node[t][time_label] = y
            for s in range(0,t):
                r= random.random()
                if r < Model.f(t,s,x,pos[s][0],c,N):
                    G.add_edge(s,t)
        if TR==True:           
            tr(G)
        return G
    
    @staticmethod
    def manhattan_distance(t_1,x_1,t_2,x_2,c=1):
        """
        Calculates Lorentzian Manhattan distance. 
        c*|t_1-t_2| + |x_1-x_2| 
        
        Input
        t_1 : time coordinate of first point
        x_1 : spatial coordinate of first point
        t_2 : time coordinate of second point
        x_2 : spatial coordinate of second point
        c=1 : scaling of time separation (speed of information)
         
        """
        if c<math.inf and c>0:
            return c* math.fabs(t_1-t_2) + math.fabs(x_1-x_2)
        return t_1-t_2
    
    @staticmethod
    def manhattan_dag(number_nodes,number_layers=1,random_on=True,c=1,d=1,TR=True,time_label='t',space_label='x'):       
        """
        Model of a DAG created using Manhattan distance.
        Total number of nodes is with number_nodes .
        Nodes are arranged in number_layers layers where each layer has the same time coordinate.
        Layers are uniformally separated with times scaled to to be between 0 to 1.
        Each x coordinate is set at uniformally random between 0 and 1.
        The edges are set if a uniformally distributed random number between 0 and d
        is greater than the Manhattan distance between the source and target node.
        
        Input
        number_nodes : number of nodes 
        number_layers=1 : Number of Layers
        random_on=True - If true place nodes randomly in spatial direction, else place uniformally distributed
        c =1 : scaling of time separation (speed of information)
        d =1 : maximum Manhattan distance allowed for connection
        TR=True : selects transitively reduced version if True
        time_label='t': Node key for time coordinate
        space_label='x': Node key for space coordinate 
        
        Return
        networkx directed graph representing Tim'source example DAG
        Nodes have keys given by time_label and space_label arguments giving time and space coordinates respectively
        
        """
        G = nx.DiGraph()
        number_per_layer = number_nodes // number_layers
        fnumber_layers = float(number_layers)
        fnumber_per_layer=float(number_per_layer)
        for target in range(number_nodes):
            G.add_node(target)
            if random_on:
                target_x = random.random()
            else: 
                target_x = (0.5+(target %number_per_layer) ) / fnumber_per_layer
            # now force integer division to split into layers
            target_t = (target//number_per_layer)/fnumber_layers
            G.node[target][space_label] = target_x
            G.node[target][time_label] = target_t
            for source in range(target):
                source_t = G.node[source][time_label]
                if target_t>source_t and d*random.random() > Model.manhattan_distance(target_t,target_x,source_t,G.node[source][space_label],c):
                    G.add_edge(source,target)
                    #print(source, '->',target)
        if TR==True:           
            tr(G)
        return G
      
     
      
    @staticmethod
    def test_dag_tim(TR=True,time_label='t',space_label='x'):        
        """
        Standard test DAG.
        This is the example used by Tim in his slides.
        Note nodes have coordinates "t" for time (a topological order) and "x" for space.
        Edges are ordered from low index nodes (and low time) to higher index nodes (and higher time).
        
        Input
        TR=True : selects transitively reduced version if True
        time_label='t': label nodes with time coordinate
        space_label='x': label nodes with space coordinate
        Return
        networkx directed graph representing Tim's example DAG
        """
    
        pos= {}
        pos[0] = (1,0)
        pos[1] = (0,1)
        pos[2] = (2,2)
        pos[3] = (0,3)
        pos[4] = (2,3)
        pos[5] = (1,4)
        G = nx.DiGraph()
        for t in range(6):
            G.add_node(t)
            G.node[t][space_label] = pos[t][0]
            G.node[t][time_label]  = pos[t][1]
        G.add_edge(0,1)    
        G.add_edge(0,2)    
        G.add_edge(1,2)    
        G.add_edge(1,3)    
        G.add_edge(2,3)    
        G.add_edge(2,4)    
        G.add_edge(3,5)    
        G.add_edge(4,5)    
        if TR==True:           
            tr(G)
        return G
