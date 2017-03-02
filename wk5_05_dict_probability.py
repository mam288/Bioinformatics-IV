"""
Solution to the Probability of Spectral Dictionary Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 5, code challenge #5
"""

import constants as co
    
def spectral_dict_size(input_vector,min_threshold,max_threshold):
    '''
    Probability of Spectral Dictionary Problem.
    Parameters
    --------
    input_vector: spectral vector (string)
    min_threshold: minimum threshold (integer)
    max_threshold:m maximum threshold (integer)
        
    Return
    --------
    Spectral Dictionary Probability (floating point number)
    '''
    vector = [0] + convert_input(input_vector)
    probability = [[0 for y in range(max_threshold+1)] for x in range(len(vector))]
    probability[0][0] = 1
    mass_values = co.masses_integers.values()
    for i in range(1,len(vector)):
        for t in range(0,max_threshold +1):
            new_probability = 0
            for mass in mass_values:
                prev_i = i-mass
                prev_t = t- vector[i]
                if prev_t <= max_threshold and prev_i >= 0 and prev_t >= 0:
                    new_probability += probability[prev_i][prev_t]
            probability[i][t] = new_probability*(1.0/len(mass_values))
    return sum(probability[-1][min_threshold:max_threshold+1])

def convert_input(input_data):
    '''
    Used with spectral_dict_size.
    '''
    split_data = input_data.split()
    vector = []
    for e in split_data:
        vector.append(int(e))
    return vector
    
###############################################################################################
if __name__ == "__main__":
    challenge_dataset = '-2 4 4 9 10 4 -9 13 2 9 -7 -2 13 -8 -6 4 1 3 9 6 -9 2 10 -1 4 11 6 -3 0 -7 12 9 -2 -3 -6 7 9 -7 10 10 1 -7 -8 13 6 -5 14 -10 -7 6 -6 6 -8 6 -2 -10 11 -5 14 -2 -10 -10 -2 -3 6 8 3 -10 -2 -10 -8 13 -4 13 4 -5 6 -8 6 -10 11 14 1 -6 5 15 13 3 4 -8 14 14 12 -2 13 -6 15 12 -6 -6 11 8 -2 15 1 5 -6 -10 -8 -1 9 9 12 -4 -9 6 -9 -6 13 15 -5 7 12 6 7 -10 -8 -10 -1 8 -6 7 -3 -5 4 10 13 4 4 -8 -5 -7 7 13 -10 6 -1 -3 5 6 -4 4 -4 7 -9 12 7 6 11 7 10 14 8 12 2 -6 9 9 4 7 2 -10 -2 -9 4 14 11 1 -2 4 4 8 -9 -9 4 -3 7 -3 8 -7 13 -8 -8 13 -7 5 14 -9 14 9 1 10 9 7 5 14 -6 -5 2 13 12 2 9 -2 -7 8 11 3 12 6 -8 9 13 -2 13 10 7 -10 -5 10 14 -3 3 4 2 12 5 -5 -2 -2 12 -10 -1 8 -9 5 -5 -5 0 -6 7 -2 -8 3 4 -6 6 10 5 -8 -9 -7 2 -6 1 -6 1 10 -8 14 10 7 -8 6 0 -9 -8 12 -4 -4 -5 8 9 1 -4 9 0 13 -5 1 1 -4 5 -8 3 1 10 -5 0 8 4 14 -7 -1 12 -9 2 5 -9 7 -2 15 -2 -5 12 -8 14 14 -2 -10 -1 5 4 -3 -1 -6 8 0 2 -8 -4 -10 6 15 -6 5 8 -5 7 3 3 1 11 12 12 4 -3 -4 -6 -10 -6 6 10 15 7 15 7 10 -8 2 5 2 9 2 9 15 11 9 14 7 10 -2 -3 9 -3 14 10 2 7 -3 -6 7 15 -7 -8 -7 6 -7 -9 12 6 -7 15 -10 -6 8 10 -7 12 -2 -7 -1 9 -3 11 2 -9 5 6 3 11 -1 4 -2 -5 -9 8 0 -6 5 -2 8 -7 3 4 13 -2 -10 -9 -10 14 -10 2 7 6 8 7 -8 1 -5 -1 10 -7 12 -5 10 -3 -8 11 14 -10 10 -6 7 -8 -8 -7 9 -2 -9 -6 10 -7 8 4 7 -1 -5 8 3 -5 -3 -4 11 3 11 -10 1 8 -8 -10 -2 9 -1 1 12 10 7 2 6 1 13 -8 14'
    min_threshold = 38
    max_threshold = 200
    print (spectral_dict_size(challenge_dataset,min_threshold,max_threshold))
