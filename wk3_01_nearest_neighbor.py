"""
Solution to the  Nearest Neighbors of a Tree Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 3, code challenge #1
"""

import copy

def neighbor_switch(input_data,node_1,node_2):
    '''
    Nearest Neighbors of a Tree Problem.
    Takes a phylogeny tree and switches two nodes adjacent to the given edge e (node1,node2).
    
    Parameters
    --------
    nodes: Phylogeny tree given as an adjacency list (string)
    node1,node2: Two internal nodes a and b specifying an edge e (integers)
    
    Return
    -------
    Nearest neighbors of the tree with respect to e.
    
    Print Output (written to results.txt): Two adjacency lists representing the neighboring trees.
    '''
    def convert_input(input_matrix):
        '''Convert input for neighbor switch'''
        split_1 = input_matrix.split('\n')
        new_dict = {}
        for row in split_1:
            row = row.split('->')
            node_1 = int(row[0])
            node_2 = int(row[1])
            if node_1 in new_dict:
                new_dict[node_1].append(node_2)
            else:
                new_dict[node_1] = [node_2]
        return new_dict
        
    def write_file(filename,content):
        with open(filename, 'a') as f:
            f.write(content + '\n')
        
    def convert_output(neighbors):
        '''Print output for neighbor_switch'''
        write_file('results_quiz.txt', str(neighbors))
        for neighbor in neighbors:
            neighbor_keys = neighbor.keys()
            neighbor_keys = [key for key in neighbor_keys]
            neighbor_keys.sort()
            for key in neighbor_keys:
                for node in neighbor[key]:
                    line = str(key) + '->' + str(node)
                    write_file('results_quiz.txt', line)
            write_file('results.txt', '\n')
            
    nodes = convert_input(input_data)
    neighbors = []
    node_1_children = nodes[node_1]
    node_2_children = nodes[node_2]
    for child_of_1 in node_1_children:
        for child_of_2 in node_2_children:
            
            # switch neighbors
            if child_of_1 != node_1 and child_of_1 != node_2 and child_of_2 != node_1 and child_of_2 != node_2:
                new_nodes = copy.deepcopy(nodes)
                
                # remove the edge between node_1 and child_of_1
                new_nodes[node_1].pop(new_nodes[node_1].index(child_of_1)) 
                new_nodes[child_of_1].pop(new_nodes[child_of_1].index(node_1))
                
                # add an edge between node_1 and child_of_2
                new_nodes[node_1].append(child_of_2); new_nodes[child_of_2].append(node_1)
                
                # remove the edge between  node_2 and child_of_2
                new_nodes[node_2].pop(new_nodes[node_2].index(child_of_2))
                new_nodes[child_of_2].pop(new_nodes[child_of_2].index(node_2))
                
                # add an edge between node_2 and child_of_1
                new_nodes[node_2].append(child_of_1); new_nodes[child_of_1].append(node_2)
                neighbors.append(new_nodes)
    convert_output(neighbors)
         
########################################################################################
if __name__ == "__main__":
    node_1 = 5
    node_2 = 4
    sample_input = '''0->4
4->0
1->4
4->1
2->5
5->2
3->5
5->3
4->5
5->4'''
    
    neighbors = neighbor_switch(sample_input,node_1,node_2)

