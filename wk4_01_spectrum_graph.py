"""
Solution to the Construct the graph of a spectrum Problem.
Molecular Evolution (Bioinformatics IV) on Coursera.org
Week 4, code challenge #1
"""

import constants as co
    
def spectrum_graph(spectrum_masses):
    '''
    Construct the graph of a spectrum.
    
    Parameters
    --------
    spectrum_masses: A space-delimited string of integers Spectrum. (string)
    
    Return
    --------
    None
    
    Print output (written to Results.txt): List of tuples representing a spectrum graph
    '''   
            
    def convert_output(spectrum_graph):
        '''Convert output for spectrum_graph.'''
        def write_file(filename,content):
            with open(filename, 'a') as f:
                f.write(content + '\n')
        
        for e in spectrum_graph:
            line = str(e[0]) + '->' + str(e[1]) + ':' + str(e[2])
            write_file('results.txt', line)
        
    spectrum_masses = convert_input(spectrum_masses)
    spectrum_graph = []
    length = len(spectrum_masses)
    mass_values = co.masses_integers.values()
    for i in range(1,length):
        for k in range(0,i+1):
            mass_diff = spectrum_masses[i] - spectrum_masses[k]
            if mass_diff in mass_values:
                spectrum_graph.append((spectrum_masses[k],spectrum_masses[i],find_key(mass_diff)))
    convert_output(spectrum_graph)

def find_key(value):
    '''Finds the amino acid associated with a given mass.'''
    for e in co.masses_integers:
        if co.masses_integers[e] == value:
            return e
    return False
    
def convert_input(input_data):
    '''Convert input for spectrum_graph.'''
    mass_list = [0]
    split_data = input_data.split()
    for e in split_data:
        mass_list.append(int(e))
    return mass_list
############################################################

if __name__ == "__main__":
    sample_input = '57 71 154 185 301 332 415 429 486'
    spectrum_graph(sample_input)
