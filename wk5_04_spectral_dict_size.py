"""
Solution to the Size of Spectral Dictionary Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 5, code challenge #4
"""

import constants as co

def spectral_dict_size(input_vector,min_threshold,max_threshold):
    '''
    Size of Spectral Dictionary Problem.
    
    Parameters
    --------
    input_vector: spectral vector (string)
    min_threshold: minimum threshold (integer)
    max_threshold:m maximum threshold (integer)
        
    Return
    --------
    Size of the dictionary consisting of all the peptides that score within the threshold limits
    against the given spectrum. (integer)
    '''
    def convert_input(input_data):
        '''
        Used with spectral_dict_size
        '''
        split_data = input_data.split()
        vector = []
        for e in split_data:
            vector.append(int(e))
        return vector
    
    vector = [0] + convert_input(input_vector)
    size = [[0 for y in range(max_threshold+1)] for x in range(len(vector))]  #matrix tracking the size of the spectral dictionary at different points along the spectrum
    size[0][0] = 1
    mass_values = co.masses_integers.values()
    for i in range(1,len(vector)):
        for t in range(0,max_threshold +1):
            new_size = 0
            for mass in mass_values:
                prev_i = i-mass
                prev_t = t- vector[i]
                if prev_t <= max_threshold and prev_i >= 0 and prev_t >= 0:
                    new_size += size[prev_i][prev_t]
            size[i][t] = new_size
    return sum(size[-1][min_threshold:max_threshold+1])

###############################################################################################
if __name__ == "__main__":
    challenge_dataset = '-5 -2 -8 9 9 11 -2 12 -1 2 7 -6 -7 1 5 6 10 -6 12 4 -7 -5 0 -6 -8 -5 7 12 7 -10 -9 -1 10 11 -7 -9 10 -1 7 13 15 13 -7 10 -10 -6 -6 8 1 9 -3 11 -7 -1 14 7 -9 4 9 -2 12 12 13 2 -9 7 -3 0 0 -8 4 -6 -9 3 -2 6 -9 -7 -1 8 9 1 9 -5 -7 1 13 -7 -2 12 -8 -3 2 -4 8 -7 12 8 5 -9 4 -6 -8 10 7 11 12 11 11 14 -4 7 11 13 12 -3 -6 -4 8 0 4 -4 9 -10 10 9 12 3 13 -1 -9 4 2 11 -6 0 9 14 7 -7 13 10 -7 8 -4 0 9 4 1 -8 14 -2 13 4 1 -5 -3 -4 15 12 9 -6 3 -3 8 -2 -4 -4 9 -4 5 6 11 -4 -10 -5 14 -5 3 -7 6 14 -10 9 -10 -1 8 13 9 -8 -9 -7 4 0 11 12 -6 -3 -10 9 -3 12 4 3 -5 -1 9 10 -5 -3 11 8 11 -9 -10 -5 -3 -2 -9 -8 -5 3 -7 10 13 -10 11 15 -9 11 8 -3 -10 -4 9 4 -2 -9 2 -2 -5 -2 -10 5 -3 -4 5 -5 3 13 -10 -6 10 15 10 -8 -9 0 -10 8 -3 9 2 0 13 -1 15 -2 10 2 5 -3 13 -7 -1 15 13 9 9 13 -1 4 -10 5 2 -9 -2 3 9 9 -5 12 -5 12 -9 -5 5 -3 -2 -1 -8 6 8 -5 15 11 5 8 4 12 -5 9 -3 9 1 6 -4 -1 9 13 -6 2 -9 14 -7 5 -8 -3 4 -4 0 14 -8 -7 0 8 12 9 14 5 -6 5 13 8 4 -7 14 1 7 5 6 4 -10 -10 -8 5 -1 9 -10 -3 6 -8 -4 11 13 4 8 9 -6 -3 -7 7 11 5 -5 11 -4 4 10 1 -4 -2 -6 -3 -10 15 8 9 8 -1 14 7 13 0 -3 8 6 15 5 13 3 4 6 -9 -7 15 6 -2 -9 -10 2 0 9 11 8 -6 9 14 1 4 3 6 9 5 14 3 0 -1 -4 5 13 1 12 12 -8 -1 -4 11 11 6 -10 -6 -8 -5 0 15 6 10 15 8 5 10 -2 -5 11 6 -9 3 10 11 9 9 3 -2 10 6 0 -2 -3 -1 7 12 -3 6 -4 13 3 14 -8 -2 -8 9 -4 -8 -9 10 0 -7 9 2 14 5 3 7 6 8'
    peptides = spectral_dict_size(challenge_dataset,35,200)
    print (peptides)
