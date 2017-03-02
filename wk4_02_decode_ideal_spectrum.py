"""
Solution to the Decoding an Ideal Spectrum Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 4, code challenge #2
"""

import constants as co
import wk4_01_spectrum_graph as sg

def decode_ideal_spectrum(input_data):
    '''
    Solve the Decoding an Ideal Spectrum Problem.
    
    Parameters 
    --------
    input_data: A space-delimited list of integers Spectrum. (string)
    
    Return
    --------
    Peptide sequence with the highest match to given spectrum (string)
    '''
    def find_ideal_spectrum_for_peptide(peptide):
        '''Calculate the ideal spectrum for a peptide.'''
        co.aa_masses = []
        ideal_spectrum = [0]
        for amino_acid in peptide:
            co.aa_masses.append(co.masses_integers[amino_acid])
        for lower_bound in range(0,len(co.aa_masses)):
            for upper_bound in range(len(co.aa_masses),lower_bound,-1):
                total = sum(co.aa_masses[lower_bound:upper_bound])
                ideal_spectrum.append(total)
        ideal_spectrum = list(set(ideal_spectrum)) # remove duplicates
        ideal_spectrum.sort()
        return ideal_spectrum
    
    def find_spectrum_graph(input_data):
        '''
        Converts a spectrum to a graph with spectrum masses as nodes and amino acids as weights. The weight is the amino acid
        with a mass equivalent to the difference between the incoming and outgoing nodes. 
        '''
        spectrum_masses = sg.convert_input(input_data)
        spectrum_graph = {}
        length = len(spectrum_masses)
        mass_values = co.masses_integers.values()
        for i in range(1,length):
            for k in range(0,i+1):
                lower = spectrum_masses[k]
                higher = spectrum_masses[i]
                mass_diff = higher - lower
                if mass_diff in mass_values:
                    amino_acid = sg.find_key(mass_diff)
                    if higher in spectrum_graph:
                        spectrum_graph[higher][lower] = amino_acid
                    else:
                        spectrum_graph[higher] = {lower:amino_acid}
        return spectrum_graph, spectrum_masses
    
    def find_paths(spectrum_graph):
        '''Finds all of the paths from source to sink in a given spectrum of masses'''
        paths = []
        masses = spectrum_graph.keys()
        masses = list(masses)
        masses.sort(reverse = True)
        sink = masses[0]                
    
        def find_paths_recursive(prev_mass, current_mass, current_path):  # m is the current mass being evaluated       
        
            if current_mass != sink:
                current_path = current_path + str(spectrum_graph[prev_mass][current_mass])
                
            if current_mass == 0:
                paths.append(current_path)
                return    
                
            if current_mass not in spectrum_graph:
                return
                
            outgoing_masses = spectrum_graph[current_mass].keys()
            
            if outgoing_masses == []:
                return 
                
            for outgoing_mass in outgoing_masses:
                find_paths_recursive(current_mass, outgoing_mass, current_path)
                
        find_paths_recursive([],sink,'')
        return paths
    
    spectrum_graph,spectrum_masses = find_spectrum_graph(input_data)
    paths = find_paths(spectrum_graph)
    for path in paths:
        ideal_spectrum = find_ideal_spectrum_for_peptide(path)
        if len(set(spectrum_masses) - set(ideal_spectrum)) == 0:
            return path
            
#######################################################################################
if __name__ == "__main__":
    
    sample_input = '57 71 154 185 301 332 415 429 486'
    
    peptide = decode_ideal_spectrum(sample_input)
    print peptide