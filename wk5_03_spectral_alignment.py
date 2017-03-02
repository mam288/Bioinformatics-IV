"""
Solution to the Spectral Alignment Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 5, code challenge #3
"""
import constants as co

def spectral_alignment(peptide,vector,max_modifications):
    '''
    Spectral Alignment Problem.
     
    Parameters
    --------
    peptide: peptide sequence (string)
    vector: spetral vetor (string)
    max_modifications: max number of modifications allowed (int)
    
    Return
    --------
    The highest scoring peptide related to the given peptide by up to max_modifications. (string)
    '''
    vector = vector.split()
    vector = [0] + [int(x) for x in vector]
    prefix_masses = calculate_mass(peptide)
    peptide_mass = prefix_masses[-1]

    # initialize the scores and backtrack matrices
    scores = [[[float('-inf') for i in range(len(vector))] for j in range(peptide_mass + 1)] for k in range(max_modifications +1)]
    backtrack = [[[float('-inf') for i in range(len(vector))] for j in range(peptide_mass + 1)] for k in range(max_modifications +1)]
    scores[0][0][0] = 0
    backtrack[0][0][0] = (0,0,0)
    
    
    for i in range(len(vector)):
        for j in prefix_masses[1:]:
            for k in range(max_modifications +1):
                prev_j_index = prefix_masses.index(j) - 1
                aa_mass = j - prefix_masses[prev_j_index]
                
                # if we are at the 0-modification level of the matrix - only consider scores from the same level (local)
                if k == 0:
                    incoming_prev_k = float('-inf')
                    local_score = max_incoming = scores[k][j-aa_mass][i-aa_mass]
                    score = local_score + vector[i]

                # if we are above the 0-modification level, consider local scores and scores from the lower level
                else:
                    if i == 0 or i-aa_mass <0 and max(scores[k-1][j-aa_mass][:i]) == float('-inf'):
                        continue
                    elif i-aa_mass < 0 and max(scores[k-1][j-aa_mass][:i]) != float('-inf'):
                        local_score = -1
                    else:
                        local_score = scores[k][j-aa_mass][i-aa_mass]

                    incoming_prev_k = scores[k-1][j-aa_mass][:i]
                    incoming_scores = incoming_prev_k + [local_score]
                    max_incoming = max(incoming_scores)
                    score = max_incoming + vector[i]

                # update the backtrack matrix
                if max_incoming == local_score:
                    backtrack[k][j][i] = (k,j-aa_mass,i-aa_mass)
                else:
                    backtrack_i = incoming_scores.index(max_incoming)
                    backtrack[k][j][i] = (k-1,j-aa_mass,backtrack_i)
                    
                scores[k][j][i] = score  # set the score for the current position

    modified_peptide = backtrack_function(backtrack,peptide)
    return modified_peptide
    
def calculate_mass(peptide):
    '''Calculates the mass of a given peptide'''
    mass_sum = 0
    prefix_masses = [0]
    for amino_acid in peptide:
        mass_sum += co.masses_integers[amino_acid]
        prefix_masses.append(mass_sum)
    return prefix_masses
    
def backtrack_function(backtrack,peptide):
    '''Backtracks to create the modified peptide which is returned as a string.'''
    modified_peptide = ''
    k = len(backtrack) - 1
    j = len(backtrack[0]) - 1
    i = len(backtrack[0][0]) - 1
    peptide_index = len(peptide) - 1
    current_amino_acid = peptide[peptide_index]

    while peptide_index >= 0:
        current_backtrack = backtrack[k][j][i]
        mass_diff = i - current_backtrack[2] - co.masses_integers[current_amino_acid]
        mass_diff_string = str(mass_diff)
        
        # if there was amino acid with a modification, calculate how much the mass differs from the original amino acid by
        if mass_diff > 0:
            mass_diff_string = '(' + '+' + mass_diff_string + ')'
        elif mass_diff <0:
            mass_diff_string = '(' + mass_diff_string + ')'
            
        # if there was no modification to the amino acid
        elif mass_diff == 0:
            mass_diff_string = ''
            
            
        modified_peptide = current_amino_acid + mass_diff_string + modified_peptide
        k = current_backtrack[0]
        j = current_backtrack[1]
        i = current_backtrack[2]
        peptide_index -= 1
        current_amino_acid = peptide[peptide_index]
    return modified_peptide


############################################################################
if __name__ == "__main__":
    code_challenge_dataset_spectrum = '-3 9 -5 3 -7 3 1 4 10 6 -1 14 -1 8 -5 2 12 -1 8 2 13 -4 15 10 14 -1 8 0 3 11 5 11 5 10 8 4 -5 8 13 -10 -2 8 5 0 -2 5 0 -8 -6 -6 12 1 -9 -4 3 4 4 9 8 9 -8 1 13 3 2 -3 3 -5 -8 -7 -9 -5 -10 -9 12 1 6 -2 6 -3 5 -3 14 -5 7 -3 -3 -1 14 3 7 1 -9 2 -9 7 7 13 7 10 10 -7 5 4 11 1 10 9 11 -10 6 5 2 -9 -9 -3 11 4 11 5 -2 3 5 10 1 -10 0 0 -1 4 -5 7 -3 15 1 -10 -10 5 5 6 10 12 2 9 5 12 11 -7 -8 13 1 -2 3 5 -3 -8 -8 15 12 2 -8 10 11 1 -7 1 -9 4 9 6 0 -6 9 15 -5 6 3 -9 3 -5 -10 -2 -6 -3 6 -7 1 9 13 -2 -10 10 -7 -2 -2 3 -10 -9 -2 12 -7 7 9 11 7 9 -5 -8 3 2 13 -1 -2 2 -8 -2 1 5 -5 10 12 -5 2 3 1 -8 -5 -9 8 -9 11 15 0 12 9 5 -9 -6 -3 4 15 10 -7 -3 -10 9 -7 12 -7 7 7 -9 -2 3 -3 -2 7 -10 -1 -9 9 10 -8 14 -4 11 15 -6 1 4 -9 2 12 7 1 -1 14 -3 1 5 -4 -4 -3 -5 -3 -8 3 -9 7 -6 15 -6 2 14 -5 -9 -3 -10 3 -6 -2 10 0 13 -10 4 8 5 0 -6 10 -4 13 11 -2 15 -7 3 13 13 -8 6 -1 5 -5 -2 7 0 -10 1 14 -7 -1 0 -4 5 7 3 -5 9 -2 1 4 4 5 -7 12 9 12 13 -1 -6 3 -8 7 -7 0 5 -3 -10 -4 -1 -8 6 15 14 -2 1 -1 0 -1 3 6 -9 -8 7 -6 6 -10 -4 9 -7 6 0 -7 10 1 4 11 -6 -9 15 7 11 -5 13 9 -5 -8 -10 -7 14 12 11 1 -4 10 6 -9 -8 13 -7 8 10 5 -9 -3 4 -3 9 1 -9 8 7 12 3 -3 5 7 6 -5 1 -3 4 -3 -6 -9 -10 10 11 12 8 13 -3 -1 2 3 -8 1 8 2 -7 -6 4 -9 10 -8 11 0 -4 2 -7 -10 -3 -4 12 -5 10 -5 -8 -1 -4 6 -8 -2 7 1 6 5 -9 12 7 -2 11 11 5 -2 8 12 12 -2 2 -9 15 -1 -8 9 5 14 -1 -4 8 -10 -4 15 -5 8 8 4 -1 5 7 2 -9 14 -10 4 -1 -5 7 -7 7 13 -5 -8 15 1 -7 12 -7 2 -5 -4 -3 6 -9 12 -7 -7 1 -9 8 -2 -5 -10 -9 -3 11 -5 -7 -4 -9 10 -6 -6 7 -7 4 -2 13 7 8 -6 7 6 11 -2 -6 -10 -1 11 -6 -8 12 4 7 5 13 -3 1 9 -4 7 -10 13 -1 -5 -6 -3 -9 -1 10 -3 3 8 -3 -8 -2 3 8 0 12 -4 -4 3 -8 13 12 -2 -7 6 8 -7 5 13 -7 13 7 13 1 13 -5 8 6 6 -6 -7 8 -3 12 10 5 9 15 1 -8 -1 11 11 -9 -4 6 -7 4 5 10 13 3 5 9 12 -2 -10 7 12 -9 1 15 -2 6 1 -5 12 -7 1 15 2 -2 15 4 1 15 9 8 1 -9 -8 -5 9 -7 -6 -6 8 -9 11 6 -9 -5 -1 1 2 -1 -6 5 10 8 0 9 -9 5 9 14 -7 7 15 -2 3 4 10 -1 11 -10 13 8 2 -8 -4 0 7 7 -6 -1 15 -2 2 -1 -7 1 13 7 12 7 3 -4 -2 3 12 6 -3 11 -1 -9 7 10 -5 -3 -6 14 3 7 6 -2 -2 3 0 -2 9 -7 8 12 3 2 10 -8 -6 6 2 13 -8 2 3 3 -3 -4 10 -10 -6 2 5 15 3 6 6 11'
    code_challenge_dataset_seq = 'STTRLIL'
    
    print spectral_alignment(code_challenge_dataset_seq,code_challenge_dataset_spectrum,3)


















