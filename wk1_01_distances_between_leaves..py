"""
Solution to Distances Between Leaves Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 1, code challenge #2
"""

def find_distances(input_data):
    '''
    Distances Between Leaves Problem. 
    
    Parameters
    --------
    input_data: Adjacency list for a phyogeny tree (string)
    
    Return
    --------
    None
    
    Print Output: A distance matrix representing the distance between all leaves on the phylogeny tree.
    '''
    def convert_input (input):
        '''Convert input for find_distances.'''
        nodes = {}
        first_split = input.split('\n')
        for node in first_split:
            split_node,distance = node.split(':')
            incoming,outgoing = split_node.split('->')
            entry = {}
            entry_2 = {}
            if incoming in nodes:
                nodes[incoming][outgoing] = int(distance)
            if outgoing in nodes:
                nodes[outgoing][incoming] = int(distance)
            if outgoing not in nodes:
                entry[outgoing] = int(distance)
                nodes[outgoing] = entry_2
            if incoming not in nodes:
                entry_2[incoming] = int(distance)
                nodes[incoming] = entry
        return nodes
    
    def find_leaves(nodes):
        '''Returns a list of leaves with no outgoing edges (leaves)'''
        leaves  = []
        for node in nodes:
            if len(nodes[node]) == 1:
                leaves.append(node)
        return leaves
               
    def traverse_graph(current, prev, visited, distance):
        '''Traverse the graph and calculate the distances between each set of leaves in the graph.'''
        visited += [current]
        if prev != '':
            distance += nodes[prev][current]
        if current in leaves and current != source:
            matrix[leaves.index(source)][leaves.index(current)] = distance
            matrix[leaves.index(current)][leaves.index(source)] = distance
            return
        neighbors = find_neighbors(current,visited)
        for neighbor in neighbors:
            traverse_graph(neighbor,current,visited,distance)
            
    def find_neighbors(origin,visited):
        '''Find unvisited neighbors'''
        neighbors = []
        for node in nodes:
            for edge in nodes[node]:
                if edge == origin and node not in visited:
                    neighbors.append(node)   
        return neighbors
        
    def format_output(matrix):
        '''Format output for find_distances'''
        for row in matrix:
            string_row = ''
            for e in row:
                string_row += str(e) + ' '
            print string_row
            
    nodes = convert_input(input_data)
    leaves = find_leaves(nodes)
    matrix = [[0]*(len(leaves)) for x in range(len(leaves))]
              
    # for each leaf traverse the graph and find the distances to all other leaves
    for source in leaves:
        traverse_graph(source,'',[],0)
        
    format_output(matrix)
    return matrix
    
##########################################################################################
if __name__ == "__main__":
    input_data = '''0->4:11
1->4:2
2->5:6
3->5:7
4->0:11
4->1:2
4->5:4
5->4:4
5->3:7
5->2:6
2->6:1
6->2:1'''

    find_distances(input_data)
