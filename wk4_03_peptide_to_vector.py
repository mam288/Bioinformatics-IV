"""
Solution to the Converting a Peptide into a Peptide Vector Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 4, code challenge #3
"""

import constants as co

def peptide_to_vector(peptide):
    '''
    Converting a Peptide into a Peptide Vector Problem.
    
    Parameters
    -----------
    peptide: Amino acid sequence (string)
    
    Return
    --------
    A vector representing the distribution of masses for the input peptide. (list of integers)
    '''
    def find_ideal_spectrum_for_peptide_prefix(peptide):
        '''Finds the ideal spectrum for a given peptide.'''
        co.aa_masses = []
        ideal_spectrum = [0]
        for amino_acid in peptide:
            co.aa_masses.append(co.masses_integers[amino_acid])
        for upper_bound in range(len(co.aa_masses),0,-1):
            total = sum(co.aa_masses[0:upper_bound])
            ideal_spectrum.append(total)
        ideal_spectrum = list(set(ideal_spectrum)) # remove duplicates
        ideal_spectrum.sort()
        return ideal_spectrum
    
    spectrum = find_ideal_spectrum_for_peptide_prefix(peptide)
    vector = [0]*max(spectrum)
    for mass in spectrum:
        vector[mass-1] = 1
    return vector
    

    
###########################################################################
if __name__ == "__main__":
    challenge_data = 'PTRCIQDSDDSSYHDDYPQFTDWVFQPDMRKLYRIYDSH'
    print (peptide_to_vector(challenge_data))
