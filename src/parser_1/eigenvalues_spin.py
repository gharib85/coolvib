import numpy as np

with open("MyM.EIG", "r") as f:
    content=f.readlines()
    fermi_level=content[0].split()[0]
    n_bands= int(content[1].split()[0])
    n_spin = int(content[1].split()[1])
    n_kpoints = int(content[1].split()[2])
    eigenvalues=np.zeros((n_bands,n_kpoints,n_spin))

    energy=[]
    counter=1.0
    line_number=2
    for k in range(0,n_kpoints):
        energy+=content[line_number].split()[1:]
        for n in range(1, int(np.ceil(n_bands*n_spin/10.0))):
            line_number+=1
            energy+=content[line_number].split()
        print energy
        for s in range (n_spin):
            for i in range(0, n_bands):
                eigenvalues[i,k,s] = energy[i+n_bands*s]
        line_number+=1
        energy=[]


    

