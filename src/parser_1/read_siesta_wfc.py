import numpy as np
import struct
import FortranFile


f= FortranFile.FortranFile("MyM.WFSX")
nk,gamma= f.readInts()
print "nk, gamma  = ", nk,gamma
nspin = int(f.readInts())
print "nspin =",nspin
nuotot = int(f.readInts())
print "nuotot =",nuotot
f.readRecord()
psi=np.zeros((nk,nspin,nuotot,nuotot))
psi = psi +0j
for iik in range (1,nk+1):  #for each k-point
    for iispin in range (1,nspin+1):  #for each spin
        f.readRecord()
        ispin = int(f.readInts())
        nwflist =int(f.readInts())
        for iw in range(1,nwflist+1):  # for each state (nwflist = total number of states)
            indwf=f.readInts()
            energy=f.readReals('d')    # we first read the energy of the state 
            read_psi = f.readReals()   # and all the orbital coefficients (real value, followed by the imaginary value
            read_psi=np.reshape(read_psi, (nwflist,2))  # reshape it 
            psi[iik-1,iispin-1,iw-1,:]=read_psi[:,0]+1j*read_psi[:,1]  # and make a row of complex numbers
print psi[55,1,70,:]


