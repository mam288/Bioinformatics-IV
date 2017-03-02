"""
Solution to the  Peptide Sequencing Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 4, code challenge #5
"""

import constants as co
    
def find_peptide_sequence(input_vector):
    '''
    Peptide Sequencing Problem.
    
    Parameters
    ---------
    input_vector: A space-delimited spectral vector 'Spectrum' (string)
    
    Return
    --------
    An amino acid string with maximum score against Spectrum (string)
    '''
    def convert_input(input_data):
        '''Convert input data for find_peptide_sequence.'''
        split_data = input_data.split()
        vector = []
        for e in split_data:
            vector.append(int(e))
        return vector
    
    def find_seq_recursive(current_mass, current_index, current_score,current_seq,track_index):    
        
        # if you have gone past the last value in the vector, the current path is not a match to the spectrum
        if current_index > vector_len - 1:
            return
            
        # if you have reached the last value in the vector 
        if current_index == vector_len-1:
            current_seq += co.reversed_masses[current_mass][0]
            current_score += vector[current_index]

            # if the current score is higher than the high score, set the current peptide as the new high-scoring peptide
            if current_score > high_score[0]:
                high_score[0] = current_score
                high_score_seq[0] = current_seq
            return
            
        # add the current value to current_score and the amino acid corresponding to the current value to the peptide
        if current_index != -1:
            current_seq += co.reversed_masses[current_mass][0]
            current_score += vector[current_index]   
     
        
        for mass in co.reversed_masses:
            if (current_index+mass) <= (vector_len - 1):
                find_seq_recursive(mass,current_index+mass,current_score,current_seq,track_index) # proceed to the next value in the vector
            else:
                return # if adding the current mass brings you beyond the end of the vector, return
                
    high_score = [float('-inf')]
    high_score_seq = ['']
    vector = convert_input(input_vector)
    vector_len = len(vector)
    find_seq_recursive(0,-1,0,'','')
    return high_score_seq[0]
    
################################################################
if __name__ == "__main__":
    challenge_data = '-3 11 -11 -8 19 19 -19 11 -18 16 17 -16 16 5 -17 3 -8 21 -11 6 -7 -14 9 25 22 8 -15 -18 -4 29 30 -1 7 -4 -19 11 17 15 -9 24 17 6 9 26 0 -4 29 -20 15 2 10 5 15 5 -4 26 2 26 12 -13 23 25 -14 19 7 -9 23 10 2 25 -7 13 -3 23 25 -18 13 2 9 26 12 13 21 15 -4 11 -15 -11 -4 -18 -13 -18 -19 7 26 18 -3 0 -6 23 29 -19 24 -4 9 21 16 7 8 -10 7 24 10 17 -14 20 25 -12 29 -12 -4 2 -9 6 -2 -4 -6 5 14 27 22 -10 27 -17 -5 -14 -20 25 8 23 -5 -1 17 -10 -14 14 17 16 -16 30 16 15 21 -7 -15 -9 -14 -16 27 24 27 19 -16 25 29 27 20 -20 16 25 -19 -9 -19 20 28 6 -2 3 2 23 -17 12 -9 -2 9 30 -15 -15 30 -12 3 23 -6 13 17 23 14 26 26 21 5 26 -18 27 -15 26 -12 12 22 0 17 -10 -3 12 18 14 1 -13 -6 -19 -1 29 18 -9 -11 -11 0 6 -17 8 5 17 15 12 -20 27 -1 -9 28 -7 19 22 -16 22 26 24 30 -12 -18 27 14 12 19 -3 20 28 -14 -8 -9 28 20 1 -9 10 22 24 23 30 2 -16 -1 20 7 -11 29 24 11 29 -9 4 -1 -6 28 11 14 17 -4 4 -4 -8 3 28 29 1 17 -13 -4 11 -3 29 11 10 -4 11 11 29 27 14 11 -10 4 -7 15 14 12 -5 13 18 -16 -8 13 -16 21 17 -1 3 16 -14 -8 5 19 12 10 3 25 -3 7 3 22 23 25 -6 7 -14 2 -4 28 6 -16 11 -3 20 -2 -2 5 6 10 5 27 30 -15 7 7 -10 23 14 21 5 -12 -13 -10 -3 23 15 -6 12 1 19 20 6 19 13 14 -11 -7 -3 -17 12 -8 -15 -4 -20 -8 30 23 -16 24 -7 4 7 18 6 6 -7 -15 -12 7 27 12 -5 19 -9 30 -10 -8 -4 -5 3 24 24 18 -12 30 -19 29 14 -7 -6 25 -8 6 14 -2 20 10 20 0 11 18 -17 22 7 12 4 29 -8 7 8 -19 9 -14 -19 -3 24 -6 -14 -18 -12 21 13 -11 27 -6 10 -19 1 -4 28 -12 13 14 -9 -10 -5 -2 25 30 -15 -16 12 -11 -3 5 -10 29 28 30 -12 16 6 26 -12 -9 -17 17 13 24 -6 18 30 -18 -6 18 6 1 -2 -8 25 -7 22 23 -12 25 10 15 19 23 -7 -16 24 -1 28 24 -18 13 -14 22 29 8 -12 -6 -7 -6 8 16 -20 24 -1 -9 7 -6 2 29 -2 6 8 10 -16 12 20 17 6 -1 -11 14 29 14 16 1 12 -11 -7 20 19 18 27 15 28 11 24 9 15 -14 7 -19 -20 -2 2 -15 14 5 0 -12 -3 10 -6 28 -18 -19 11 -18 5 -15 -1 -20 27 3 5 8 5 -15 -19 -12 14 -16 6 -20 23 8 1 16 -7 24 16 -8 22 18 -14 6 10 3 -14 -4 11 29 -10 9 -15 6 22 -12 14 -8 12 -4 10 23 12 26 -8 0 7 -13 10 -12 29 22 -15 5 15 27 -17 25 -3 -14 20 5 -14 6 -11 1 13 9 10 -11 -3 -14 30 -19 30 -20 6 -4 -8 -16 -20 30 11 -18 -5 -14 4 21 28 29 -5 16 7 4 30 12 15 12 -17 21 -9 0 1 10 -7 9 -12 -11 6 18 12 -3 21 -8 9 -17 13 14 23 -12 29 -15 9 -16 -18 -4 5 17 14 7 -4 -16 17 30 6 11 -14 -10 -12 4 -20 9 21 6 -3 -11 17 -4 -3 -15 -16 -17 -6 21 26 -14 17 -14 -16 -3 26 13 -1 -19 -17 -10 28 10 10 15 -12 15 18 -19 7 -4 -7 -13 9 7'
    print find_peptide_sequence(challenge_data)

