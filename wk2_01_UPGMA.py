"""
Implementaton of UPGMA
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 2, code challenge #1
"""

import wk1_02_limb_length as ll
import wk1_03_additive_phylogeny as ap

def upgma(matrix):
    '''
    UPGMA
    
    Parameters
    --------
    matrix: Space separated m x m distance matrix. (string)
    
    Return
    ---------
    None
    
    Print Output: Adjacency list for the phylogeny tree created using UPGMA
    '''
    
    def set_ages(matrix):
        '''Set the age for each node in the graph to 0'''
        ages = {}
        for i in range(len(matrix) + 1):
            ages[i] = 0.0
        return ages
        
    def merge_clusters(nodes, clusters,cluster_a,cluster_b,m):
        '''Merge cluster a and cluster b'''
        clusters[m] = {}
        for cluster in clusters:
            clusters[m][cluster] = ()
            dist_a,size_a = clusters[cluster][cluster_a][0], clusters[cluster][cluster_a][1]
            dist_b, size_b= clusters[cluster][cluster_b][0],clusters[cluster][cluster_b][1]
            size_total = size_a + size_b
            avg_dist = (dist_a*size_a + dist_b*size_b)/size_total
            clusters[cluster][m] = clusters[m][cluster] = (avg_dist,size_total)
            if cluster == cluster_a:
                edge_dist = dist_b/2
                nodes = add_edge(nodes,cluster_a,cluster_b,edge_dist,m)
            if m == cluster: 
                clusters[cluster][m] = (0.0,0)
            clusters[cluster].pop(cluster_a); clusters[cluster].pop(cluster_b)
        m += 1
        clusters.pop(cluster_a), clusters.pop(cluster_b)
        return clusters,nodes,m
    
    def add_edge(nodes,node_a,node_b,edge_dist,m):
        '''Add an edge to the graph'''
    
        def create_new(node1,node2,dist):
            '''Create a new node and add an edge'''
            nodes[node1] = {node2:dist}

        def add_to_existing(node1,node2,dist):
            '''Add an edge using an existing node'''
            nodes[node1][node2] = dist

        def type_specific_addition(node1,node2):
            '''Determine whether to add an edge to an existing node or create a new node'''
            if node2 in nodes and len(nodes[node2]) != 1:
                dist = edge_dist - ages[node2]
            elif node1 in nodes and len(nodes[node1]) != 1:
                dist = edge_dist - ages[node1]
            else:
                dist = edge_dist
            if node1 not in nodes:
                create_new(node1,node2,dist)  
            else:
                add_to_existing(node1,node2,dist)
        
        # add edges between node_a to m and m to node_b
        type_specific_addition(node_a,m); type_specific_addition(m,node_a)
        type_specific_addition(node_b,m); type_specific_addition(m,node_b)
        return nodes    
                
    clusters = create_clusters(matrix,False)
    matrix = ll.convert_input(matrix)
    m = len(matrix)
    ages = set_ages(matrix)
    nodes = {}

    # while more than 1 cluster exists find the closest nodes and merge into a cluster
    while len(clusters) > 1:
        cluster_a, cluster_b,dist = find_closest_nodes(clusters)
        ages[m] = dist/2
        clusters, nodes,m = merge_clusters(nodes,clusters,cluster_a,cluster_b,m)
    
    ap.convert_output(nodes)
            
def find_closest_nodes(clusters):
    '''Find the closest clusters in the graph'''
    min_dist = float('inf')
    min_cluster = ()
    for cluster in clusters:
        for entry in clusters[cluster]:
            current_dist = clusters[cluster][entry][0]
            if current_dist < min_dist and cluster != entry:
                min_dist = current_dist
                min_cluster = (cluster,entry,min_dist)
    return min_cluster

def create_clusters(matrix,is_neighbor_joining):
    '''Create a dictionary of clusters using the given matrix'''
    matrix = ll.convert_input(matrix)
    clusters = {}
    for i in range(len(matrix)):
        for j in range(len(matrix)):
            if i not in clusters:
                if is_neighbor_joining == False:
                    clusters[i] = {j:(float(matrix[i][j]),1)}   # set the size of the cluster to 1 and value to matrix value
                else:
                    clusters[i] = {j:(float(matrix[i][j]))}
            else: 
                if is_neighbor_joining == False:
                    clusters[i][j] = (float(matrix[i][j]),1)
                else:
                    clusters[i][j] = (float(matrix[i][j]))
    return clusters
                
##############################################################################
if __name__ == "__main__":
    sample_input = '''0	20	17	11
    20	0	20	13
    17	20	0	10
    11	13	10	0'''
    
    nodes = upgma(sample_input)


