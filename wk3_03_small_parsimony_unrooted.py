"""
Solution to Small Parsimony in an Unrooted Tree Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 3, code challenge #3
"""

import networkx as nx
import numpy as np
import wk3_02_small_parsimony as sp

def small_parsimony_unrooted(input_data):
    '''
    Small Parsimony in an Unrooted Tree Problem..
    
    Parameters
    --------
    input_data: Adjacency list for an unrooted binary tree (string)
        
    Return
    --------- 
    The minimum parsimony score of this tree 
    
    Print output: Adjacency list of the small parsimony tree
    '''

    def calc_min(graph):    
        '''Calculate the minimum parsimony for the graph'''
        root = max(graph.nodes())
        total = 0
        for nuc in graph.node[root]['score']:
            total += min(nuc)
        return total,root
        
    def calculate_score(graph,node):
        '''Calculate the score and update the backtrack attribute for the given node'''
        parent_score = graph.node[node]['score']
        children = graph.successors(node)
        for i in range(len(parent_score)):
            min_child_scores = {}
            min_indices = {}
            score_total = {}
            child_scores = {}

            # get the score for each child and store it in child_scores
            for child in children:
                child_score = graph.node[child]['score'][i]   # get the score for each child of the current node
                min_child_score = min(child_score) 
                child_scores[child] = child_score
                min_child_scores[child] = min_child_score
            
            # calcuate the total score for all children - including necessary mismatch penalties 
            for n in range(4):  
                for child in children:
                    
                    # find the indices where the minimum score occurs 
                    child_min_indicies = graph.node[child]['backtrack'][i] = np.where(child_scores[child] == min_child_scores[child])
                    min_indices[child] = child_min_indicies

                    # add mismatch penalties if appropriate 
                    if n in min_indices[child][0]:
                        mismatch = 0
                    else:
                        mismatch = 1
                        
                    score_total[child] = min_child_scores[child] + mismatch

                total_all_children = sum(score_total.values()) # total score for all chidren
                graph.node[node]['score'][i][n] =  total_all_children

        graph.node[node]['tag'] = 1
        return graph
        
    def assemble_seqs(graph,root):
        '''Use the graph to create the minimum parsimony sequences at each node'''
        nucs = ['A','C','G','T']
        nodes_list = graph.nodes()
        nodes_list.sort()
        seq_len = len(graph.node[root]['score'])
        for i in range(seq_len):
            min_score = min(graph.node[root]['score'][i])
            score_list = graph.node[root]['score'][i]
            index = np.where(score_list == min_score)[0][0]
            graph.node[root]['seq'] += nucs[index]
            graph.node[root]['backtrack'][i] = [index]
        nodes_list = nodes_list[:-1]
        for node in nodes_list[::-1]:
            if graph.node[node]['seq'] != '':
                continue
            if graph.predecessors(node) == []:
                continue
            parent = graph.predecessors(node)[0]
            for i in range(seq_len):
                
                # find the nucleotides with the lowest scores at the parent node for the current position
                # this will help determine which nucletide to assign to the current position
                parent_chars_set = np.unique(graph.node[parent]['backtrack'][i]) 
                
                # find the nucleotides with the lowest scores at the current position (the assigned nucleotide has to be chosed from this group)
                child_chars_set = np.unique(graph.node[node]['backtrack'][i]) 
                
                # find out which nodes are on the lowest score lists for both the parent and current nodes
                possible_assignments = np.intersect1d(parent_chars_set,child_chars_set) 
                possible_assignments = sorted(possible_assignments)
                
                # if there are no nucletides common to both lists, pick a nucleotide that is on the lowest score list for the child node
                if possible_assignments == []:  
                    index = sorted(child_chars_set)[0]

                # else, pick one of the nucleotides that is on both lowest score lists
                else:   
                    index = possible_assignments[-1]
                graph.node[node]['seq'] += nucs[index]
                graph.node[node]['backtrack'][i] = [index]
        return graph
        
    def find_ripe_node(graph):   # determine if there are any ripe nodes in the tree
        for node in graph.nodes():
            if graph.node[node]['tag'] == 1:
                continue
            children = graph.successors(node)
            if children == []:
                continue
            for child in children:
                if graph.node[child]['tag'] == 0:
                    continue
            return node
        return False
     
    def add_root(graph,seq_len):
        '''Add a root to the unrooted graph between the highest numbered and second highest numbered nodes.'''
        largest_node_num = max(graph.nodes()) # highest numbered node
        root = largest_node_num + 1
        graph.add_node(root,tag = 0,seq = '',score = np.full((seq_len,4),float('inf')),backtrack = ['b' for x in range(seq_len)])
        
        # find the secon highest numbered node
        outgoing_nodes_from_largest = graph.successors(largest_node_num)
        second_largest_node = outgoing_nodes_from_largest[::-1][0]

        # remove the existing edge and add new outgoing edges from the root
        graph.remove_edge(largest_node_num,second_largest_node)
        graph.add_edges_from([(root,largest_node_num),(root,second_largest_node)])
        return graph, root, (largest_node_num,second_largest_node)
        
    def remove_root(graph,root,removed_edge):
        '''Remove the root from the tree and add back the edge that was removed'''
        graph.remove_node(root)
        graph.add_edge(removed_edge[0],removed_edge[1])
        return graph
        
    def remove_backwards_edges(graph):
        '''Remove all edges in the graph that do not go in the root->leaf direction'''
        for incoming,outgoing in graph.edges():
            if incoming < outgoing:
                graph.remove_edge(incoming,outgoing)
        return graph
            
    def create_graph(input_data):
        '''Create and initialize the graph, including the nodes, scores, seqs, and backtrack attributes for each node.'''
        
        def create_score(graph,node,seq):
            '''Use the sequence to calculate the score for the given node'''
            seq_len = len(seq)
            nucs = ['A','C','G','T']
            for i in range(seq_len):
                index = nucs.index(seq[i])
                graph.node[node]['score'][i][index] = 0
            return graph
            
        def add_digit(graph,new_node,seq_len):
            '''Add a node to the graph given an integer'''
            graph.add_node(new_node,tag = 0,seq = '',score = np.full((seq_len,4),float('inf')),backtrack = ['b' for x in range(seq_len)])
            return graph
        
        def add_seq(graph,seq,node):
            '''Add a node to the graph given a sequence'''
            graph.add_node(node,tag = 1,seq = seq,score =np.full((len(seq),4),float('inf')),backtrack = ['b' for x in range(len(seq))])
            graph = create_score(graph,node,seq)
            return graph
            
        graph = nx.DiGraph()
        data = input_data.split('\n')
        leaf_num = 0
        seq_len = len(data[0].split('->')[0])
        seq_dict = {}
        for e in data:
            entry = e.split('->')
            incoming = entry[0]

            # if the incoming node is an integer
            if entry[0].isdigit():
                incoming = int(entry[0])
                graph = add_digit(graph,incoming,seq_len)
                
            # if the incoming node is a sequence
            else:
                if entry[0] not in seq_dict.keys():
                    incoming = leaf_num
                    graph = add_seq(graph,entry[0],incoming)
                    seq_dict[entry[0]] = incoming
                    leaf_num += 1
                else:
                    incoming = seq_dict[entry[0]]

            outgoing = entry[1]

            # if the outgoig node is an integer            
            if entry[1].isdigit():
                outgoing = int(entry[1])
                graph = add_digit(graph,outgoing,seq_len)
                
            # if the outgoing node is a sequence
            else:
                if entry[1] not in seq_dict.keys():
                    outgoing = leaf_num
                    graph = add_seq(graph,entry[1],outgoing)
                    seq_dict[entry[1]] = outgoing
                    leaf_num += 1
                else:
                    outgoing = seq_dict[entry[1]]

            graph.add_edge(incoming,outgoing)
        return graph,seq_len
        
    # if the input data is a string, create the graph
    if type(input_data) == str:
        graph,seq_len = create_graph(input_data) 
        graph = remove_backwards_edges(graph)
        
    # if the input data is already a graph, set the graph variable and seq_len
    else:
        graph = input_data; seq_len = len(graph.node[0]['backtrack'])
        
    graph,root,removed_edge = add_root(graph,seq_len)
    ripe_node = find_ripe_node(graph)
    
    # while there are still unvisited nodes calculate the score for 'ripe' nodes
    while ripe_node != False:   
        graph = calculate_score(graph,ripe_node) # calculate the scores for the ripe node
        ripe_node = find_ripe_node(graph) 
        
    min_parsimony,root = calc_min(graph)
    
    # use the calculated scores to assemble the sequences with the lowest parsimony score at each node
    graph = assemble_seqs(graph,root)  
    
    graph = remove_root(graph,root,removed_edge)
    sp.convert_output(graph)
    print min_parsimony
    return graph, min_parsimony
                
################################################################
if __name__ == "__main__":
    sample_input = '''TCGGCCAA->4
4->TCGGCCAA
CCTGGCTG->4
4->CCTGGCTG
CACAGGAT->5
5->CACAGGAT
TGAGTACC->5
5->TGAGTACC
4->5
5->4'''
    small_parsimony_unrooted(sample_input)



