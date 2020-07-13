# Siblinarity-based antichain partitioning of DAGs

**Siblinarity-based antichain partitioning of DAGs**

Code for antichains written by Vaiva Vasiliauskaite v.vasiliauskaite16@imperial.ac.uk and Tim Evans t.evans@imperial.ac.uk .

Vasiliauskaite, V., Evans, T.S. Making communities show respect for order. Appl Netw Sci **5**, 15 (2020). https://doi.org/10.1007/s41109-020-00255-5
---

## Source Files

1. node_matrix_greedy_antichain_partition.py
	This is the main file, which contains function matrix_node_recursive_antichain_partition. The function takes into its argument a graph,
	whose partition will be found and returned. The method is analogous to Louvain community detection algorithm.
	Step one is to place each node into its own community. After that, for each node in the network, its second 
	nearest neighbours (predecessors' successors, SNNs) are found and it is placed into the antichain with its SNN,
	for which the Quality function increase is the largest (given it is larger than 0). After that, the graph is 
	coarse-grained: each antichain is a supernode and edges between antichains are edges between supernodes. After
	the new graph is built, the algorithm comes back to step one. The loop is continued until the quality function
	no longer increases. 
		
2. model_antichain_dag.py
	This is a model of a spatial DAG to investigate the quality of the partitioning algorithm. In the model, two nodes
	u,v connect with an edge (u,v) with a probability p, proportional to dist(u,v). pos(u) = x_u,y_u; the distance 
	is calculated as Manhattan distance. If y_u>y_v, edge (u,v) never occurs.
	
3. utilities_antichains.py
	The file contains all supplementary functions used in the antichain_partition algorithm.
	is_antichain - Tests whether a list of nodes in a DAG is an antichain.
  tr - finds transitively reduced version of a DAG
  is_weakly_connected - tests whether a list of source nodes in a graph have a path in either direction between a list of target nodes.
  
4. Quality_matrix
  Quality measures for use in antichains. Similarity matrix implementation. Definition is that used for weighted graph.

  \begin{equation}Q = \sum_{u \in partition1} \sum_{v \in partition2}( S_ij- k_i k_j/W )\end{equation}
  where W = total strength of edges in the graph (($sum_{i,j}S_ij)/2$),
 $ S_{ij}$ - i,j^th entry in the similarity matrix. For instance, $A.A^T$ is successors-based similarity;
          $A^T.A$ is predecessors-based similarity.   
	  

