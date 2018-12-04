# -*- coding: utf-8 -*-
"""
Created on Mon Oct  1 10:20:54 2018
@author: Vaiva
***********IMPORTANT NOTE***************
Only implemented for directed networks without cycles. For cyclic networks, use node_greedy_antichain_partition
We can only use numpy implementation on DAGs. The coarse-grained network is likely to be
cyclic, so it is better to use matrices only in the first step, where 
the heavy linear algebra lifting is done. After this, we use the usual weakly_connected_graph antichain check.
"""

from quality_matrix import *
from utilities_antichain import *
import datetime
import random 
from model_antichain_dag import Model
import itertools
import math
import numpy as np
import networkx as nx
from collections import defaultdict
import scipy.sparse as sparse

def coarse_grain(G,node_to_partition_label, partition_label_to_nodes, weight_attribute='weight',
                 time_label='t',space_label='x'):
    ''' Coarse Grain
    
    The new graph H has the partitions of G as the nodes in H.
    An edges from partition1 to partition2 in H is present if
    there is an edge from a node in G in partition1 of G to 
    a node in G in partition2.  The total weight of the edge 
    from partition1 to partition2 in H will be the sum of all 
    the weights of all such edges in G 
    from nodes in partition1 to nodes in partition2.
    If unweighted, weights are assumed to be 1. 
    If time_label or space_label are set, these are assumed to be numerical 
    values (e.g. coordinates) and nodes in the new graph get the average value from
    the partition of nodes they represent in the old graph.
    
    
    Input
    ----
    G - networkx graph
    node_to_partition_label - dictionary from G node key to its partition label
    partition_label_to_nodes - dictionary from partition label to the set of nodes in G in that partition
    weight_attribute='weight' - attribute on edge containing edge weight data
    time_label='t': Node key for time coordinate (used as y/vertical coordinate)
    space_label='x': Node key for space coordinate (used as x/horizontal coordinate)
    
    Return
    ------
    H - coarse grained graph, nodes are the partitions labels, weights are under eight_attribute of edges
    
    '''
    H = nx.DiGraph()
    H.add_nodes_from(list(partition_label_to_nodes.keys()) )
    for partition in partition_label_to_nodes.keys():
        nodes_in_partition = partition_label_to_nodes[partition]
        number_in_partition = len(partition_label_to_nodes[partition])
        if time_label!=None:
            average_time = sum([get_node_attribute_value(G,n,node_attribute=time_label) for n in partition_label_to_nodes[partition] ]) / number_in_partition 
        H.node[partition][time_label]=average_time
        if space_label!=None:
            average_space = sum([get_node_attribute_value(G,n,node_attribute=space_label) for n in partition_label_to_nodes[partition] ]) / number_in_partition 
        H.node[partition][space_label]=average_space
                
        
    for partition1, partition2 in itertools.combinations(partition_label_to_nodes.keys(),2):
        w= sum( [ get_edge_weight(G,node1,node2, weight_attribute)  for node1, node2  in itertools.product(partition_label_to_nodes[partition1], partition_label_to_nodes[partition2] ) ] )
        if w>0:
            H.add_edge(partition1,partition2,weight_attribute=w)                
        w= sum( [ get_edge_weight(G,node2,node1, weight_attribute)  for node1, node2  in itertools.product(partition_label_to_nodes[partition1], partition_label_to_nodes[partition2]) ] )
        if w>0:
            H.add_edge(partition2,partition1,weight_attribute=w)
    return H    
            

def similarity_matrix_sparse(DAG, similarity = "intersection",neighbours = "successors"):
    """
    Function to produce a sparse similarity matrix based on neighbourhoods of nodes in DAG.
    
    Input
    -----
    DAG - networkx directed acyclic graph
    similarity - type of similarity of sets. Currently only implemented for the size of intersection
    neighbours - type of neighbours to consider in the similarity. Can be either successors or predecessors
    
    Return
    -----
    A - scipy sparse symmetric similarity matrix where entry A[i,j] represents similarity between nodes of indices i, j 
    nodedict - dictionary of node names and their indices in the similarity matrix
    """
    
    nodes = list(DAG.nodes())
    nodedict = {}
    for i in range(len(nodes)):
        nodedict[nodes[i]] = i
    nodelist = list(nodedict.keys())
        
    A = (nx.adjacency_matrix(DAG,nodelist))
    
    if similarity == "intersection" and neighbours == "successors":
        A = A.dot(A.transpose())
        A.setdiag(0)
        return A,nodedict
    
    elif similarity == "intersection" and neighbours == "predecessors":
        A=A.transpose().dot(A)
        A.setdiag(0)
        return A,nodedict


def similarity_matrix(DAG, similarity = "intersection",neighbours = "successors"):
    """
    Function to produce a similarity matrix based on neighbourhoods of nodes in DAG.
    
    Input
    -----
    DAG - networkx directed acyclic graph
    similarity - type of similarity of sets. Currently only implemented for the size of intersection
    neighbours - type of neighbours to consider in the similarity. Can be either successors or predecessors
    
    Return
    -----
    A - symmetric similarity matrix where entry A[i,j] represents similarity between nodes of indices i, j 
    nodedict - dictionary of node names and their indices in the similarity matrix
    """
    
    nodes = list(DAG.nodes())
    nodedict = {}
    for i in range(len(nodes)):
        nodedict[nodes[i]] = i
    nodelist = list(nodedict.keys())
        
    A = (nx.adjacency_matrix(DAG,nodelist)).todense()
    
    if similarity == "intersection" and neighbours == "successors":
        A = A.dot(A.transpose())
        np.fill_diagonal(A,0)
        return A,nodedict
    
    elif similarity == "intersection" and neighbours == "predecessors":
        A=(A.transpose()).dot(A)
        np.fill_diagonal(A,0)
        return A,nodedict
    
    

def has_path_matrix(DAG,nodedict,cutoff = 350):
    
    nodes = list(nodedict.keys())
    A = (nx.adjacency_matrix(DAG,nodes)).todense()
    A_sum = np.copy(A)
    if nx.is_directed_acyclic_graph(DAG)==True:
        L_max = len(nx.dag_longest_path(DAG))
    else:
        L_max = cutoff
    current_length = 1
    #for i in range(2,L_max):
    #    A_sum += (np.linalg.matrix_power(A,i)>0).astype(np.int_)
    while current_length <= L_max:
        A_sum= np.dot(A_sum,A)
        current_length +=1
        
    return (A_sum >0).astype(np.int8)#,nodedict


def find_paths_sparse(A,length_max=10):
    '''
    Tim's numpy path implementation updated by Vaiva to sparse matrices
    Scipy sparse matrix implementation to find all paths
    Given adjacency matrix A will find all the paths between all vertices
    
    Input
    A: numpy square adjacency matrix,can be weighted
    Return
    #tuple current_length,path_length,path_bool where
    #current_length = one more than the maximum length found.  
                   If equals length_max then may have terminated because reachd maximum requested length
    #non_zero_entries = number of non-zero entries in P=(A)^current_length
    #path_length = matrix of longest paths lengths 
                   path_length[target,source]= longest path from source to target
    path_bool = matrix of booleans indicating if path exists. Paths from vertex to slef (length zero) gives True on diagonal
                   path_bool[target,source]= True (False) if path from source to target
    '''
    
    # Assume vertices start connected to selves only path path of length zero.
    m,n= np.shape(A)
    path_bool=sparse.eye(m,n,dtype = bool)
    path_bool=path_bool.tolil()
    path_length= sparse.lil_matrix((m, n), dtype = "int32")
    current_length=1
    # P(n+1)=P(n)*A 
    P=A.copy()
    non_zero_entries=P.count_nonzero()

    while non_zero_entries>0 and current_length<length_max:
        non_zero_entries=P.nonzero()
        path_bool[non_zero_entries[0],non_zero_entries[1]]=True
        path_length[non_zero_entries]=current_length
        P=P.dot(A)#np.dot(P,A)
        current_length += 1
        non_zero_entries=P.count_nonzero()
    #return current_length,path_bool
    return path_bool


def is_weakly_connected_matrix(path_matrix, nodedict, source_nodes,target_nodes):
    """
    Function to check whether nodes in the source_nodes are not weakly connected to
    nodes in the target_nodes. 
    
    Input
    -----
    path_matrix - 1/0 matrix where entry P[i,j] = 1 if nodes with indices i,j are weakly connected
    nodedict - dictionary where keys are node names and values are their corresponding indices in the path matrix
    source_nodes - list of nodes
    target_nodes - list of nodes
    
    Return
    ------
    True - if nodes in source_nodes and target_nodes form a weakly_connected subgraph
    False - if not
    """
    source_nodes_id , target_nodes_id = [nodedict[s] for s in source_nodes],[nodedict[t] for t in target_nodes]
    
    for s,t in itertools.product(source_nodes_id,target_nodes_id):
        if path_matrix[s,t] == 1 or path_matrix[t,s] == 1:
            return True
    return False


def node_matrix_greedy_antichain_partition(G,level,
                            random_on=False,seed=None, max_number_sweeps=None,
                            backwards_forwards_on=True,forwards_backwards_on=False,
                            Q_check_on=True, weight_attribute = "weight"):    
    """
    In this implementation we iterate over nodes in the graph moving individual nodes until no changes occur .
    
    We start with each node in its own partition.
    In one sweep we look at each partition ac in turn. 
    We find all the backwards-forwards neighbours of the nodes in partition ac 
    and collect all their partition labels, excluding the current partition ac.
    For each of these we find if we can increase the quality function by merging 
    current partition ac with one of its bf-neighbouring partitions. If we do we then
    do the merge removing the current partition ac.
    We continue the sweep looking at remaining partitions and trying to merge them.  Note that
    later partitions will see some previous partitions already merged.  The option to randomise the order
    that we visit the partitions will lead to different results.
    
    After each sweep, if at least one partition was merged then we sweep through again.
    We only stop when no more merges are found on one sweep, or if the number of sweeps exceed
    the maximum requested.
    
    Note that we do NOT create an explicit hybrid graph. We only use a weakly connected check
    of the  
    
    Input
    -----
    G - networkx graph. 
    random_on=False - if true will shuffle the order in which partitions are examined
    seed=None - used as seed if shuffle is on.  If None then time is used as seed
    max_number_sweeps=None - this is the maximum number of sweeps to consider. If less than 1 or None, then uses number of nodes.
    backwards_forwards_on=True  - find possible new partitions by making a backwards step then a forwards step from node being considered for a move
    forwards_backwards_on=False - find possible new partitions by making a forwards step then a backwards step from node being considered for a move
    Q_check_on=False - check to see if change in Q is correct by printing out total value and changes
    weight_attribute - edge attribute of weight. if None, unweighted quality functions are used. Note, weight must be integer
    
    Return
    ------
    tuple node_to_partition_label, partition_label_to_nodes 
    where
    node_to_partition_label is a dictionary from node key to its partition label
    partition_label_to_nodes is a dictionary from partition label to the set of nodes in that partition
   
    """ 
    if not (forwards_backwards_on or backwards_forwards_on):
        raise  ValueError("At least one of forwards_backwards_on or backwards_forwards_on parameters must be True")

    
    
    if backwards_forwards_on == True:
        #use in-degree quality with backwards-forwards step (predecessors)
        adj_matrix, nodedict = similarity_matrix(G, similarity ="intersection",neighbours="predecessors")
    elif forwards_backwards_on== True:
        #use out-degree quality with forwards-backwards step (successors)
        
        adj_matrix, nodedict = similarity_matrix(G, similarity ="intersection",neighbours="successors")
    Q = Quality_matrix(nodedict,adj_matrix)
    Q_method = Q.delta_strength_quality_unnormalised
    Q_total = Q.total_strength_quality_unnormalised

    number_of_nodes = G.number_of_nodes()
    
    if max_number_sweeps == None or max_number_sweeps <1:
        max_number_sweeps = number_of_nodes
    if random_on:
        random.seed(seed)
    # set up partition node dictionaries
    # These play the role of induced graphs
    # first place each node into its own partition            
    node_to_partition_label = {}
    partition_label_to_nodes = {}
    next_partition_label=0
    
    for n in G:
        node_to_partition_label[n]=next_partition_label
        partition_label_to_nodes[next_partition_label]=set([n])
        next_partition_label+=1
    moved=True
    number_sweeps = 0

    if Q_check_on:
        Q_total_current= Q_total(partition_label_to_nodes.values())
        print (' Initial Q value is '+str(Q_total_current))
    
    if level==0:
        L_max = len(nx.dag_longest_path(G))
        path_matrix = find_paths_sparse(nx.adjacency_matrix(G,list(nodedict.keys())),L_max)
        #path_matrix = find_paths(nx.adjacency_matrix(G,list(nodedict.keys())).todense(),L_max)#path_matrix = has_path_matrix(G,nodedict)#(adj_matrix,l_max)
    count= 0
    
    while moved == True and number_sweeps < max_number_sweeps:
        # Start of one sweep through all current partitions.
        # Check every partition ac in turn, doing the best merge you can for each partition ac.
        # Note the partition ac under study will be removed in the merge.
        # That means the list of partition labels will be altered and therefore 
        # this list can not be used as a list of partition labels in the ac loop.
        # For that reason we need a deep copy of the current list of partition labels
        # which we get from the keys of the partition label to node set dictionary
        # Conveniently we can shuffle the list used in teh ac loop if we want to randomise the greedy
        # algorithm
        number_sweeps+=1
        number_moves=0
        node_list = list(G.nodes())
        if random_on:
            random.shuffle(node_list)    
        for n in node_list :
            count+= 1
            #print(count)
            # check to see if node n should be moved
            moved = False
            # ac is the partition containing node n
            ac=node_to_partition_label[n]
            # now find the contribution from Q that comes if we move n into its own partition
            partition_ac_no_n = set(partition_label_to_nodes[ac])
            partition_ac_no_n.discard(n)
            if len(partition_ac_no_n)>0:
                delta_Q_remove_n = Q_method([n],partition_ac_no_n )
            else:
                delta_Q_remove_n = 0
            
            # now find the neighbouring partitions via backwards-forward step
            bf_nearest_neighbours_all = set()
            if backwards_forwards_on:
                for p in G.predecessors(n):
                    bf_nearest_neighbours_all.update(G.successors(p))
            if forwards_backwards_on:
                for p in G.successors(n):
                    bf_nearest_neighbours_all.update(G.predecessors(p))
            bf_nearest_neighbour_partitions =  set(node_to_partition_label[bf_nn] for bf_nn in bf_nearest_neighbours_all)
            # remove current partition ac from neighbours
            try:
                bf_nearest_neighbour_partitions.remove(ac)
            except KeyError: # ac is not in set, must have no neighbours or no b-f nearest neighbours
                pass    
            
            # dictionary from partition label to delta quality value, 
            # so delta_Q_nn[ac_nn] is change in quality if node n was to join partition ac_nn
            delta_Q_nn = {}
            # Loop round bf nearest neighbour partitions ac_nn.
            # Check each ac_nn partition to make sure it is not weakly connected to partition ac
            # then calculate the modularity change if partitions ac and ac_nn are merged
            for ac_nn in bf_nearest_neighbour_partitions:
                if level ==0:
                    if not  is_weakly_connected_matrix(path_matrix,nodedict,[n],partition_label_to_nodes[ac_nn] ) :
                        delta_Q_nn[ac_nn] = Q_method([n], partition_label_to_nodes[ac_nn])
                        #print("level",level,"using path matrix", "no of nodes =" , G.number_of_nodes())
                else:
                    if not is_weakly_connected(G,[n],partition_label_to_nodes[ac_nn] ):
                        delta_Q_nn[ac_nn] = Q_method([n], partition_label_to_nodes[ac_nn])
                        #print("level",level,"using networkx", "no of nodes =" , G.number_of_nodes())
            if len(delta_Q_nn)> 0:
                # Note nice use of operator.itemgetter to get key with largest value
                # https://stackoverflow.com/questions/268272/getting-key-with-maximum-value-in-dictionary
                ac_max = max(delta_Q_nn, key=delta_Q_nn.get)
                if delta_Q_nn[ac_max] > delta_Q_remove_n and ac_max != ac:
                    # now merge partition ac into partition ac_max 
                    number_moves+=1
                    node_to_partition_label[n]=ac_max
                    partition_label_to_nodes[ac_max].add(n)
                    partition_label_to_nodes[ac].remove(n)
                    if len(partition_label_to_nodes[ac])==0:      # no more elements in this partition
                        partition_label_to_nodes.pop(ac,None) # remove ac from this dictionary

                    if Q_check_on:
                        dQ=delta_Q_nn[ac_max] - delta_Q_remove_n
                        print('Merging node '+str(n)+' into partition '+str(ac_max)+' for increase in Q of '+str(dQ))
                        Q_total_old=Q_total_current
                        Q_total_current= Q_total(partition_label_to_nodes.values())
                        _error=Q_total_current-Q_total_old-dQ
                        if math.fabs(_error)<1e-6:                            
                            print (' Now Q='+str(Q_total_current)+', consistentcy OK ' )
                        else:
                            print (' Now Q='+str(Q_total_current)+', *** BAD consistentcy '+str(_error) )
                    moved = True
                    
                elif delta_Q_remove_n <0:
                    print('We could allow node '+str(n)+' to form a new partition for an increase in Q of '+str(- delta_Q_remove_n))
                    number_moves+=1
                    node_to_partition_label[n]=next_partition_label
                    partition_label_to_nodes[next_partition_label]=set([n])
                    next_partition_label+=1
                    
                
        # keeping looping through all partitions until can not merge any more 
    # keep doing new  sweeps as long as something changed
    return node_to_partition_label, partition_label_to_nodes    



def matrix_node_recursive_antichain_partition(G,time_label='t',space_label='x',
                            random_on=False,seed=None, max_number_sweeps=None,
                            backwards_forwards_on=True,forwards_backwards_on=False,
                            Q_check_on=True,
                            plot_on=False,
                            filenameroot=None, extlist=['pdf'],
                            ScreenOn=False):
    '''
    Use , **kwargs in func defn and call with kawargs the dictionary for named arguments
    used for the partition
    '''

    result_list=list()

    _matrix_node_recursive_antichain_partition_step(G,
                            time_label=time_label,space_label=space_label,
                            level=0,
                            result_list=result_list,
                            random_on=random_on,
                            seed=seed,
                            max_number_sweeps=max_number_sweeps,
                            backwards_forwards_on=backwards_forwards_on,
                            forwards_backwards_on=forwards_backwards_on,
                            Q_check_on=Q_check_on,
                            plot_on=plot_on,
                            filenameroot=filenameroot,extlist=extlist,
                            ScreenOn=ScreenOn)

    return result_list

    
def _matrix_node_recursive_antichain_partition_step(G,time_label='t',space_label='x',
                                             level=0,result_list=None,
                                             random_on=False,
                                             seed=None, max_number_sweeps=None,
                            backwards_forwards_on=True,forwards_backwards_on=False,
                            Q_check_on=True, 
                            plot_on=False,
                            filenameroot=None, extlist=['pdf'],
                            ScreenOn=False):    
    '''Internal routine to perform recursive version of node greedy
    '''
    print ('--- level'+str(level))
    result_list.append(None)
    node_to_partition_label, partition_label_to_nodes = node_matrix_greedy_antichain_partition(G,
                            random_on=random_on,
                            seed=seed, 
                            max_number_sweeps=max_number_sweeps,
                            backwards_forwards_on=backwards_forwards_on,
                            forwards_backwards_on=forwards_backwards_on,
                            Q_check_on=Q_check_on,
                            level = level)
    if len(partition_label_to_nodes.keys()) == G.number_of_nodes():
        return 
    # optional plot
    if plot_on:
        #ScreenOnValue=True
        node_labels_on_value=True
        cluster_labels_on_value=True
        plot_dag(G,time_label=time_label,space_label=space_label,
#             filenameroot=None,extlist=None,
#             messageString='Plot',ScreenOn=ScreenOnValue,
             filenameroot=filenameroot+'_l{0:d}'.format(level),
             extlist=extlist,
             messageString='Plot',
             ScreenOn=ScreenOn,
             node_labels_on =node_labels_on_value,
             cluster_dict=node_to_partition_label,
             cluster_labels_on=cluster_labels_on_value,
             cluster_to_nodes=partition_label_to_nodes)    

    new_G = coarse_grain(G,node_to_partition_label, partition_label_to_nodes)
    _matrix_node_recursive_antichain_partition_step(new_G,
                            time_label=time_label,space_label=space_label,
                            level=level+1,
                            result_list=result_list,
                            random_on=random_on,
                            seed=seed, 
                            max_number_sweeps=max_number_sweeps,
                            backwards_forwards_on=backwards_forwards_on,
                            forwards_backwards_on=forwards_backwards_on,
                            Q_check_on=Q_check_on, 
                            plot_on=plot_on,
                            filenameroot=filenameroot,
                            extlist=extlist,
                            ScreenOn=ScreenOn
                            )  
    result_list[level] = {'level':level, 'n_to_p':node_to_partition_label, 'p_to_n':partition_label_to_nodes }
    return 

if __name__ == "__main__":  
    G = Model.spatial_antichain_dag(100,10)
    result_list = matrix_node_recursive_antichain_partition(G)
    node_partition = {}
    
    for n in G.nodes():
    
        p_at_level= result_list[0]["n_to_p"][n]
        for i in range(1,len(result_list)-1):
            p_at_level = result_list[i]["n_to_p"][p_at_level]
        node_partition[n]=p_at_level
        
    partition = {}
    partition = defaultdict(list)
    for key,val in node_partition.items():
        partition[val].append(key)
