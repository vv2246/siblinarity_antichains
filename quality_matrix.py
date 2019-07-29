# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 16:38:51 2018
@author: Vaiva
"""

import itertools
class Quality_matrix:
    '''
    Quality measures for use in antichains. Similarity matrix implementation.
    '''
    # class variables
    def __init__(self,node_id_dict,similarity_matrix,Lambda):
        '''Initial Quality Measures
        '''
        self.similarity_matrix = similarity_matrix
        self.node_id_dict = node_id_dict
        self.strength = similarity_matrix.sum(axis = 0)
        self.strength = {n:self.strength[0,node_id_dict[n]] for n in node_id_dict.keys()}#sum over rows?
        self.total_weight = similarity_matrix.sum()#/2
        self.Lambda = Lambda
        
    def delta_strength_quality_unnormalised(self,partition1,partition2):
        '''Using in-strength null model calculate the change in unnormalised quality if two partitions are combined.
        
        Definition is that used for weighted graph.
        
        Q = \sum_{u \in partition1} \sum_{v \in partition2}
            ( S_ij
             - k_i*k_j/W )
        where W = total strength of edges in the graph ((sum_{i,j}S_ij)/2),
        S_{ij} - i,j^th entry in the similarity matrix. For instance, A.A^T is successors-based similarity;
                A^T.A is predecessors-based similarity.        
        
        Note this is not normalised.  
        
        Note no test for connectedness of nodes in partitions.
        
        Note both partitions must be non-empty otherwise TypeError raised.
        
        Note no test to see if partitions are sets or if they share common elements.
        
        Input
        partition1 - iterable list or set of the nodes in first partition
        partition2 - iterable list or set of the nodes in second partition
        
        Return
        Contribution of the quality Q from the all pairs of nodes with one from partition1, second from partition2 
        '''
        S = 0
        for node1 in partition1:
            for node2 in partition2:
                #if node1!= node2:
                S+=self.similarity_matrix[self.node_id_dict[node1],self.node_id_dict[node2]]- self.Lambda*self.strength[node1]*self.strength[node2]/self.total_weight 
        
        return S
        #return sum([self.similarity_matrix[self.node_id_dict[node1],self.node_id_dict[node2]]
        #           - self.Lambda*self.strength[node1]*self.strength[node2]/self.total_weight 
        #           for node1, node2  in itertools.product(partition1, partition2) ] )
        
    
        
    def total_strength_quality_unnormalised(self,partitions):
        '''Calculate the total unnormalised quality using strength null model
        
        Definition is that used for weighted graph.
        
        Q = \sum_{u \in partition1} \sum_{v \in partition2}
            ( S_ij
             - k_i*k_j/W )
        where W = total strength of edges in the graph ((sum_{i,j}S_ij)/2),
        S_{ij} - i,j^th entry in the similarity matrix. For instance, A.A^T is successors-based similarity;
                A^T.A is predecessors-based similarity.        
        
        Note this is not normalised.  
        
        Note no test for connectedness of nodes in partitions.
        
        Input
        -----
        partition - iterable list of sets of nodes in partitions so that
                    partition[p] is an set (or any iterable list) of the nodes in partition p
        
        Return
        ------
        Total value of the quality Q from the all pairs of nodes 
        '''

        S= 0
        for p in partitions:
            for node1 in p:
                for node2 in p:#itertools.combinations(p,2):
                    #print(self.similarity_matrix[self.node_id_dict[node1],self.node_id_dict[node2]], "-",self.strength[node1],"x",self.strength[node2],"by",self.total_weight )
                
                    S+=(self.similarity_matrix[self.node_id_dict[node1],self.node_id_dict[node2]]  - self.Lambda*self.strength[node1]*self.strength[node2]/(self.total_weight))

        return S
        #return sum( [ 
        #                sum( [  self.similarity_matrix[self.node_id_dict[node1],self.node_id_dict[node2]]  
        #                - self.Lambda*self.strength[node1]*self.strength[node2]/(self.total_weight)
        #                      for node1, node2 in itertools.combinations(p,2)] ) 
        #            for p in partitions ] )
    
    
    
