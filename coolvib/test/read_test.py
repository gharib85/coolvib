
from sys import exit

if __name__=="__main__":
    pass
else:
    exit()


from fortran_file_class import FortranFile
f=FortranFile("MyM.WFSX")
n_kpoints,gamma=f.readInts("i")
print n_kpoints, gamma
nspin=int(f.readInts("i"))
print nspin
nbasis=int(f.readInts("i"))
print nbasis
f.readRecord()
import numpy as np
#coeff_real=np.zeros([nbasis,n_kpoints,nspin])
#coeff_complex=np.zeros([nbasis,n_kpoints,nspin])
coeff=np.zeros([nbasis,nbasis,n_kpoints,nspin],dtype=np.complex)

#for i in range(n_kpoints):
for i in range(1):
    for s in range(nspin):
        f.readRecord() 
        ispin=f.readInts("i")
        nwflist=f.readInts("i")
        print ispin
        print nwflist
        for j in range(int(nwflist)):
#        for i in range(5):
            indwf = f.readInts("i")
            #print indwf
            energy = f.readReals("d")
            #print energy
            tmp=None
            tmp=f.readReals("f")
            tmp=tmp.reshape(-1,2)
            print tmp
            #coeff_real= tmp[:,0]
            #coeff_complex=tmp[:,1]
            coeff[:,j,i,s]=tmp[:,0] +1.j*tmp[:,1]
            print np.dot(tmp[:,0]+1j*tmp[:,1],tmp[:,0]-1j*tmp[:,1])
            #print coeff[:,j,i,s]
            #print np.dot(coeff[:,j,i,s], coeff[:,j,i,s].conj())



