"""
Solution to the Small Parsimony Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 3, code challenge #2
"""

import networkx as nx
import numpy as np

def small_parsimony(input_data):
    '''
    Small Parsimony Problem.
    
    Parameters
    --------
    input_data: Adjacency list for an unrooted binary tree whose n leaves are labeled by DNA strings and
     whose internal nodes are labeled by integers. (string)
        
    Return
    --------- 
    The minimum parsimony score of this tree (floating point number)
    
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
        child_1,child_2 = graph.successors(node)

        for i in range(len(parent_score)):
            
            # get the score for each child of the current node
            child_1_score = graph.node[child_1]['score'][i]   
            child_2_score = graph.node[child_2]['score'][i]
            
            # find minimum score among the children
            min_child_1 = min(child_1_score); min_child_2 = min(child_2_score)  
            
            # find the indices where the minimum score occurs for each child
            child_1_min_indicies = graph.node[child_1]['backtrack'][i] = np.where(child_1_score == min_child_1)
            child_2_min_indicies = graph.node[child_2]['backtrack'][i] = np.where(child_2_score == min_child_2)
            
            # add mismatch penalities when applicable
            for n in range(4):  
                if n in child_1_min_indicies[0]: # if the nucletides match, do not add a mismatch penalty
                    mismatch_child_1 = 0
                else:    # if the nucletides don't matrch, add a mismatch penalty of 1
                    mismatch_child_1 = 1
                if n in child_2_min_indicies[0]:
                    mismatch_child_2 = 0
                else:
                    mismatch_child_2 = 1
                    
                # total the scores for each child
                total_child_1 = min_child_1 + mismatch_child_1
                total_child_2 = min_child_2 + mismatch_child_2
                
                # add the totals for both children and set that as the score for that position
                graph.node[node]['score'][i][n] = total_child_1 + total_child_2 

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
                
                # if there are no nucleotides common to both lists pick a nucleotide from the chile's low score list
                if possible_assignments == []:  
                    index = sorted(child_chars_set)[0]

                # else, pick one of the nodes that is on both lowest score lists
                else:   
                    index = possible_assignments[-1] # assign a nucleotide from the lowest score list for the current node
                    
                graph.node[node]['seq'] += nucs[index]
                graph.node[node]['backtrack'][i] = [index]

        return graph
        
    def find_ripe_node(graph):   
        '''Find a ripe node in the tree. Return False if there are no ripe nodes'''
        for node in graph.nodes():
            if graph.node[node]['tag'] == 1:
                continue
            children = graph.successors(node)
            for child in children:
                if graph.node[child]['tag'] == 0:
                    continue
            return node
        return False
            
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
            
        graph = nx.DiGraph()
        data = input_data.split('\n')
        leaf_num = 0  # number the leaves starting with 0
        seq_len = len(data[0].split('->')[1])
        
        for e in data:
            entry = e.split('->')
            incoming = int(entry[0])
            graph.add_node(incoming,tag = 0,seq = '',score = np.full((seq_len,4),float('inf')),backtrack = ['b' for x in range(seq_len)])
            
             # if both the incoming and outgoing nodes are integers
            if entry[0].isdigit() and entry[1].isdigit(): 
                outgoing = int(entry[1])
                graph.add_node(outgoing,tag = 0,seq='',score = np.full((seq_len,4),float('inf')),backtrack = ['b' for x in range(seq_len)])
            
            # if the outgoing node is a sequence
            else:
                outgoing = leaf_num
                seq_len = len(entry[1])
                graph.add_node(outgoing,tag = 1,seq = entry[1],score =np.full((seq_len,4),float('inf')),backtrack = ['b' for x in range(seq_len)])
                graph = create_score(graph,outgoing,entry[1])
                leaf_num += 1
            graph.add_edge(incoming,outgoing)
        return graph,seq_len

    graph,seq_len = create_graph(input_data)
    ripe_node = find_ripe_node(graph)
    
    while ripe_node != False:   # while there are still unvisited nodes calculate the score for 'ripe' nodes
        graph = calculate_score(graph,ripe_node) # calculate the scores for the ripe node
        ripe_node = find_ripe_node(graph)
        
    min_parsimony,root = calc_min(graph)
    graph = assemble_seqs(graph,root)  # use the calculated scores to assemble the sequences with the lowest parsimony score at each node
    convert_output(graph)
    return min_parsimony

def convert_output(graph):
        '''Convert output for small_parsimony'''
        def hamming_dist (str1, str2):
            '''Finds the number of mismatches between two strings of the same length'''
            mismatches = 0
            for i in range (0,len(str1)):
                if str1[i] != str2[i]:
                    mismatches += 1
            return mismatches
            
        for node in graph.nodes():
            children = graph.successors(node)
            for child in children:
                line1 = str(graph.node[node]['seq']) + '->' + str(graph.node[child]['seq']) + ':' + str(hamming_dist(graph.node[node]['seq'],graph.node[child]['seq']))
                line2 = str(graph.node[child]['seq']) + '->' + str(graph.node[node]['seq']) + ':' + str(hamming_dist(graph.node[node]['seq'],graph.node[child]['seq']))
                print(line1);print(line2)
                


################################################################
if __name__ == "__main__":
    sample_input = '''4->CAAATCCC
4->ATTGCGAC
5->CTGCGCTG
5->ATGGACGA
6->4
6->5'''
    
    print small_parsimony(sample_input)

