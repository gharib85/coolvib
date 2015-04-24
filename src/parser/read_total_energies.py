import numpy as np


"""
    Reads the total energy from the siesta output file 
    The way the utility is writen now the file will be called 
    input.out
"""
filename = "input"
with open("%s.out" % filename, "r") as f:
    for line in f.readlines():
        if "siesta:         Total"  in line:
            total_energy = line.split()[3]

print total_energy

#return total_energy
