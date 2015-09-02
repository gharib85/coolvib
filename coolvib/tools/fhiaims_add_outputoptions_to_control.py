"""
This script takes an FHI-AIMS output generated with 
keyword 'output k_point_list' and prints the necessary 
input parameters for the control.in file that enable writing 
of Hamiltonian, Overlap and eigenvectors.

Execute this script the following way::

    fhiaims_add_outputoptions_to_control.py <OUTPUT_filename>

OUTPUT_filename has to be an FHI-AIMS OUTPUT file generated with control.in keyword
'output k_point_list'
"""
#!/usr/bin/python

from sys import argv, exit
import numpy as np

if __name__=="__main__":

    if len(argv)<2:
        print 'Execute this script the following way:'
        print 'fhiaims_add_outputoptions_to_control.py <OUTPUT_filename>'
        print 'OUTPUT_filename has to be an FHI-AIMS OUTPUT file generated with keyword '
        print "'output k_point_list '"
    else:
        print 'Reading FHI-AIMS Output file'

        filename = argv[1]
        lines = open(filename,'r').readlines()
        kpoints = []
        for l,line in enumerate(lines):
            #if 'k_point_list' in line:
                #pass
            #else:
                #raise RuntimeError('OUTPUT file has to have been generated with \
                        #keyword: output k_point_list.')

            if 'K-points in task  ' in line:
                nkpts = int(line.split()[-1])
                n_line = l+1
                while 'K-points in task  ' in lines[n_line]:
                    nkpts += int(lines[n_line].split()[-1])
                    n_line += 1
                for nk in range(nkpts):
                    kpoint = []
                    kpoint.append(float(lines[n_line].split()[4]))
                    kpoint.append(float(lines[n_line].split()[5]))
                    kpoint.append(float(lines[n_line].split()[6]))
                    kpoint.append(float(lines[n_line].split()[9]))
                    kpoints.append(kpoint)
                    n_line += 1
                break

        kpoints = np.array(kpoints)

        print kpoints
        print ''
        n_kpts = len(kpoints)
        n_bands = int(n_kpts/2)
        
        output = '#add this to your control.in file'
        output += 'output k_point_list \n'
        
        for n in range(n_bands):
            output += 'output band  '
            for i in range(2):
                for xyz in range(3):
                    output += ' {0:8.4f} '.format(kpoints[n*2+i,xyz])
            output += ' 2 \n'

        output += 'output eigenvectors \n' 
        output += 'output overlap_matrix \n' 
        output += 'output hamiltonian_matrix \n' 

        print output

else:
    pass
