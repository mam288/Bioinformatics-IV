"""
Implementaton of  the NeighborJoining Algorithm.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 2, code challenge #2
"""
import wk2_01_UPGMA as upgma
import copy

def neighbor_joining_outer(matrix, n,m):
    '''
    NeighborJoining
    
    Parameters
    --------
    matrix: n x n distance matrix (string)
    n: matrix size (integer)
    m: keeps track of internal node number (integer)
    
    Return
    --------
    None
    
    Print Output: Adjacency list for a phylogeny tree create from the given distance matrix.
    '''
    def neighbor_joining_inner(matrix, n,m):
        '''Inner neighbor_joining function'''
        
        def construct_neighbor_joining_matrix(matrix):  
            '''Construct n-j matrix using total distances.'''
            neighbor_joining_matrix = copy.deepcopy(matrix)
            total_dist = [0 for x in range(max(matrix.keys())+1)]
            for i in matrix.keys():
                total_dist[i] = sum(matrix[i].values())
            for i in matrix.keys():
                for j in matrix.keys():
                    if i == j:
                        continue
                    neighbor_joining_matrix[i][j] = (n-1)*matrix[i][j] - total_dist[i] - total_dist[j]
            return neighbor_joining_matrix,total_dist
            
        def find_closest_nodes(matrix):  
            '''Find neighboring nodes.'''
            min_dist = float('inf')
            min_nodes = ()
            for i in matrix.keys():
                for j in matrix.keys():
                    current_dist = matrix[i][j]
                    if current_dist < min_dist and i != j:
                        min_dist = current_dist
                        min_nodes = (i,j,min_dist)
            return min_nodes
            
        def shrink (matrix,node_a,node_b):   
            '''Combine rows and shrink matrix to create D*'''
            new_index = max(matrix.keys())+1
            new_matrix = copy.deepcopy(matrix)
            for k in matrix.keys() + [max(matrix.keys())+1]:
                if k in matrix:
                    new_dist = (matrix[k][node_a] + matrix[k][node_b] - matrix[node_a][node_b])/2
                else:
                    new_dist = 0
                if new_index in new_matrix:
                    new_matrix[new_index][k] = new_dist
                else:
                    new_matrix[new_index] = {k:new_dist}
                new_matrix[k][new_index] = new_dist
            new_matrix.pop(node_a)
            new_matrix.pop(node_b)
            for entry in new_matrix:
                new_matrix[entry].pop(node_a); new_matrix[entry].pop(node_b)
            return new_matrix,new_index
        
        def add_limbs(nodes,node_a,node_b,new_index,limb_length_a,limb_length_b):
            '''Add limbs to the tree.'''
            
            def create_new(node1,node2,limb_length):
                '''Create a new node and add a limb.'''
                nodes[node1] = {node2:limb_length}

            def add_to_existing(node1,node2,limb_length):
                '''Add a limb to an existing node.'''
                nodes[node1][node2] = limb_length

            # add limb fromo node a to new node
            if node_a not in nodes:
                create_new(node_a,new_index,limb_length_a)
            else:
                add_to_existing(node_a,new_index,limb_length_a)
                
            # add limb from node b to new node
            if node_b not in nodes:
                create_new(node_b,new_index,limb_length_b)
            else:
                add_to_existing(node_b,new_index,limb_length_b)
                
            # add limbs from new node to node a and and node b
            if new_index not in nodes:
                create_new(new_index,node_a,limb_length_a)
                add_to_existing(new_index,node_b,limb_length_b)
            else:
                add_to_existing(new_index,node_a,limb_length_a)
                add_to_existing(new_index,node_b,limb_length_b)
            
        if n == 1:      #base case - return the tree with the initial limb created
            nodes = {}
            key1, key2 = matrix.keys()
            value = matrix[key1][key2]
            nodes[key1] = {key2:value}
            nodes[key2] = {key1:value}
            return nodes
        n_j_matrix,total_dist = construct_neighbor_joining_matrix(matrix) # construct n-j matrix 
        node_a,node_b,min_dist = find_closest_nodes(n_j_matrix)  # find neighboring nodes
        dist = (total_dist[node_a]-total_dist[node_b])/(n-1)
        limb_length_a = (matrix[node_a][node_b]+ dist)/2; limb_length_b = (matrix[node_a][node_b]- dist)/2  #calculate the limb_length
        matrix,new_index = shrink(matrix,node_a,node_b)   # shrink into D*
        nodes = neighbor_joining_inner(matrix,n-1,m+2)  # repeat neighbor_joining with the new matrix
        add_limbs(nodes,node_a,node_b,new_index,limb_length_a,limb_length_b) # go back and add the limbs to the tree using the distances calculated in each round of neighbor_joining
        return nodes
    matrix = upgma.create_clusters(matrix,True)
    upgma.convert_output(neighbor_joining_inner(matrix,n,m))

            
##############################################################################
if __name__ == "__main__":
    sample_input = '''0	23	27	20
    23	0	30	28
    27	30	0	30
    20	28	30	0'''
    
    neighbor_joining_outer(sample_input,3,4)







