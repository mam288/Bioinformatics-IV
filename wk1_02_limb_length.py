"""
Solution to the Limb Length Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 1, code challenge #2
"""



def find_limb(matrix,j):
    '''
    Limb Length Problem.
    
    Parameters
    --------
    matrix: Space-separated additive distance matrix D (string or list of lists)
    j: Row number (integer)
        
    Return
    --------
    The limb length of the leaf in Tree(D) corresponding to row j of this distance (integer)
    '''
    if type(matrix) == str:
        matrix = convert_input(input_matrix)
    if len(matrix) == 2:  # if the matrix is 2x2, return the distance between the two leaves
        return matrix[0][1]
    row_j = matrix[j]
    min_len = float('inf')
    for i in range(len(row_j)):
        for k in range(len(row_j)):
            if i != j and k != j and i != k: 
                current_len = (row_j[i] + row_j[k] - matrix[i][k] )/2
                if current_len < min_len:
                    min_len = current_len
    return min_len

def convert_input(input_matrix):
    '''Convert input for find_limb'''
    split_1 = input_matrix.split('\n')
    split_2 = []
    for row in split_1:
        row = row.split()
        row = [int(x) for x in row]
        split_2.append(row)
    return split_2
    
###########################################################33
if __name__ == "__main__":
    input_matrix = '''0	13	21	22
    13	0	12	13
    21	12	0	13
    22	13	13	0'''
    print (find_limb(input_matrix,1))