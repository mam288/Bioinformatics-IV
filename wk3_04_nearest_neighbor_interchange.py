"""
Implementation of the nearest neighbor interchange heuristic for the Large Parsimony Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 3, code challenge #4
"""

import copy
import wk3_02_small_parsimony as sp
import wk3_03_small_parsimony_unrooted as spu
import numpy as np
    
def nearest_neighbor_interchange(input_data,n):
    '''
    Use the nearest neighbor interchange heuristic algorithm to solve the large parsimony problem.
       
    Parameters
    --------
    input_data: Adjacency list for an unrooted binary tree with n leaves (string)
    n: number of leaves in the given tree (integer)
     
    Return
    --------
    Neighboring tree with the lowest parsimony score. (NetworkX DiGraph)
    
    Print Output: The parsimony score and unrooted labeled tree obtained after every step of the nearest neighbor interchange heuristic.
    '''
    
    def neighbor_switch(graph,node_1,node_2):
        '''Finds a set of the nearest neighbors for a phylogeny tree relative to a given edge (node1,node2).'''
        neighbors = []

        # find the outgoing neighbors for node_1 and node_2
        node_1_children = graph.successors(node_1)
        node_2_children = graph.successors(node_2)
        
        for child_of_1 in node_1_children:
            for child_of_2 in node_2_children:
                
                # switch neighbors
                if child_of_1 != node_1 and child_of_1 != node_2 and child_of_2 != node_1 and child_of_2 != node_2:
                        new_graph = copy.deepcopy(graph)
                        
                        # add new edges from node_1 and node_2 to their new neighbors
                        new_graph.add_edge(node_2,child_of_1)
                        new_graph.add_edge(node_1,child_of_2)
                        
                        # remove old edges
                        new_graph.remove_edge(node_1,child_of_1)
                        new_graph.remove_edge(node_2,child_of_2)
                        
                        # append the new graph to the list of neighbors
                        neighbors.append(new_graph)
        return neighbors
    
    def reset(graph):
        '''Re-initialize the score, backtrak, seq and tag values for the internal nodes of the graph'''
        seq_len = len(graph.node[0]['backtrack'])
        internal_nodes = [x for x in graph.nodes_iter() if graph.out_degree(x)!=0]
        for node in internal_nodes:
            graph.node[node]['score'] = np.full((seq_len,4),float('inf'))
            graph.node[node]['backtrack'] = ['b' for x in range(seq_len)]
            graph.node[node]['seq'] = ''
            graph.node[node]['tag'] = 0
        return graph
           
    score = float('inf')
    graph, min_parsimony = spu.small_parsimony_unrooted(input_data)  
    new_score = min_parsimony
    new_tree = graph
    num_nodes = len(graph.nodes())
    
    while new_score < score:
        score = new_score
        current_tree = new_tree
        
        # search the neighboring trees for a new minimum parsimony score
        for i in range(n,num_nodes):
            for neighbor_node in current_tree.nodes():
                if i in current_tree.neighbors(neighbor_node):
                    
                    # find the neighboring trees
                    neighbor_trees = neighbor_switch(current_tree,i,neighbor_node)
                    
                    # for each tree on the neighboring trees list, calculate the minimum parsimony score
                    for neighbor_tree in neighbor_trees:
                        neighbor_tree = reset(neighbor_tree)
                        new_graph, new_min_parsimony = spu.small_parsimony_unrooted(neighbor_tree)
                        
                        # set new_score to the new min parsimony score if the new min parsimony score is smaller than the current score
                        if new_min_parsimony < new_score:
                            new_score = new_min_parsimony
                            new_tree = new_graph
        sp.convert_output(graph)
        print new_score
    return new_tree
     
#################################################################################################
if __name__ == "__main__":
    sample_input = '''GCAGGGTA->5
TTTACGCG->5
CGACCTGA->6
GATTCCAC->6
5->TTTACGCG
5->GCAGGGTA
5->7
TCCGTAGT->7
7->5
7->6
7->TCCGTAGT
6->GATTCCAC
6->CGACCTGA
6->7'''
    print (nearest_neighbor_interchange(sample_input,5))