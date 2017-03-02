"""
Implementaton of AdditivePhylogeny to solve the Distance-Based Phylogeny Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 1, code challenge #3
"""

import wk1_02_limb_length as ll
import copy

def additive_phylogeny_recursive(matrix,n,m):
    """
    Implement AdditivePhylogeny to solve the Distance-Based Phylogeny Problem.
    
    Parameters
    --------
    matrix: n x n distance matrix (string)
    n: matrix size (integer)
    m: keeps track of internal node number (integer)
        
    Return
    --------
    None 
    
    Print Output:  A weighted adjacency list for the simple tree fitting this matrix.    
    """
    def check_for_dist(nodes,node_x,dist):
        '''
        Returns (True,outgoing node) if the graph has an edge of length distance leading to or from node_x. Returns (False,outgoing node)
        if there is an edge but it is not of the correct distance. Returns (False, None) if node_x is not in nodes at all. 
        '''
        if node_x not in nodes:
            return (False,None)
        else:
            for edge in nodes[node_x]:
                if (nodes[node_x][edge] == dist):
                    return (True,edge)
                else:
                    return (False,edge)
        return (False,None)

    def i_to_k_path(nodes,current,final,path):
        '''Returns a list of the nodes that are along the path from node i to node k'''
        visited = []
        final_path = []
        def rec_func(nodes,current,final,path):
            path = path + [current]
            visited.append(current)
            neighbor_nodes = nodes[current].keys()
            if current == final:
                final_path.extend(path) # Mutate external scoped list with final path value
                return
                
            unvisited_neighbor_nodes = set(neighbor_nodes) - set(visited)
            if len(unvisited_neighbor_nodes) == 0:
                return # Reached a leaf we don't care about
            for unvisited_neighbor_node in list(unvisited_neighbor_nodes):
                rec_func(nodes,int(unvisited_neighbor_node),final,path)
    
        rec_func(nodes,current,final,path)
        return final_path
    
    def find_i_k(matrix,n):
        '''Find indices for nodes i,k so node n can be inserted along this path.'''
        bald_matrix = matrix
        for i in range(len(bald_matrix)):
            for k in range (len(bald_matrix)):
                if (bald_matrix[i][k] == bald_matrix[n][k] + bald_matrix[n][i]) and (k != n) and (i != n):
                    return (i,k)
                    
    def add_to_graph(nodes,n):
        """
        Add limbs to the graph.
        Has access to current i, k, x, matrix vals within scope of recursion
        """
        
        def add_to_node(next_node):
            """Function to add limb to existing node"""
            edge_dist = ll.find_limb(matrix,n)
            nodes[next_node][n] = edge_dist
            nodes[n] = {next_node:edge_dist}
            return nodes
            
        def split_graph_between(cur_node,next_node,dist1,dist2):
            """Function to split graph and add limb to leaf from new node"""
            edge_dist = ll.find_limb(matrix,n)
            
            # Disconnect existing edge
            nodes[cur_node].pop(next_node)
            nodes[next_node].pop(cur_node)
            
            # Insert new node in between
            nodes[cur_node][m[0]] = dist1
            nodes[next_node][m[0]] = dist2
            nodes[m[0]] = {cur_node:dist1, next_node:dist2}
                  
            # Insert new leaf connected to m
            nodes[m[0]][n] = edge_dist
            nodes[n] = {m[0]:edge_dist}
            m[0] = m[0] + 1
            return nodes
            
        i_k_path = i_to_k_path(nodes,i,k,[])   #find the path frpm i to k
        total_dist = 0
        
        # add limb along i-k path
        for cur_indx in range(len(i_k_path)):
            
            # calculate the distance to the next node along the path
            cur_node = i_k_path[cur_indx]
            next_node = i_k_path[cur_indx+1]
            dist_between = nodes[cur_node][next_node]
            total_dist += dist_between
            
            if total_dist == x:
                return add_to_node(next_node)  # add limb to to an existing node
            elif total_dist > x:
                dist1 = x - (total_dist - dist_between)
                dist2 = total_dist - x
                return split_graph_between(cur_node,next_node,dist1,dist2) # create a new node to add the limb
        return nodes
        
    if n == 1:
        nodes = {}
        nodes[1] = {0:matrix[0][1]}
        nodes[0] = {1:matrix[0][1]}
        return nodes
    limb_length = ll.find_limb(matrix,n)
    new_matrix = copy.deepcopy(matrix)
    for j in range(n):
        
        # subtract limb length from all distances to and from n in the matrix
        new_matrix[j][n] = new_matrix[j][n] - limb_length
        new_matrix[n][j] = new_matrix[j][n]   

    (i,k) = find_i_k(new_matrix,n)
    x = new_matrix[i][n]

    # remove row n and column n
    new_matrix.pop(n)   
    for row in new_matrix:
        row.pop(n)
        
    nodes = additive_phylogeny_recursive(new_matrix,n-1,m)
    nodes = add_to_graph(nodes,n)
    return nodes
                
def convert_output(nodes):
    '''Convert output for additive_phylogeny'''
    nodes_keys = nodes.keys()
    nodes_keys = [int(key) for key in nodes_keys]
    nodes_keys.sort()
    for key in nodes_keys:
        for node in nodes[key]:
            line = str(key) + '->' + str(node) + ':' + str(format(nodes[key][node], '.3f'))
            print line

###########################################################33

if __name__ == "__main__":
    input_matrix = '''0	13	21	22
    13	0	12	13
    21	12	0	13
    22	13	13	0'''
        
    matrix = ll.convert_input(input_matrix)    
    tree = additive_phylogeny_recursive(matrix,3,[4])    
    convert_output(tree)
