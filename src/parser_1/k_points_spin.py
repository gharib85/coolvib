import numpy as np

with open("MyM.KP", "r") as f:
    content=f.readlines()
    n_kpoints=int(content[0].split()[0])
    kpoints_weights=np.zeros((n_kpoints))
    line_number=1
    while line_number <=len(content)-1:
        kpoints_weights[line_number-1]=content[line_number].split()[4]
        line_number+=1

