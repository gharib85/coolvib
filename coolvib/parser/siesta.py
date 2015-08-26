"""
Parsing routines for Siesta

Copyright Mikhail Askerka and Reinhard Maurer, Yale University, 03/17/2015
"""

import numpy as np
import struct
from scipy.io import FortranFile
import siesta_mod
from vibcooling.constants import bohr,ryd

def siesta_read_eigenvalues(filename):
    """
    This routine reads the number of K points, the eigenvalues and the \
    and the number of states
    Output is given in units of eV.

    Input: the name of the eigenvalue file
    Output: Fermi_level, Eigenvalues i

    The Eigenvalues is a 3D numpy array:
    1) K point
    2) Spin
    3) State

    written by Mikhail Askerka, Yale University, 03/17/2015
    """

    with open("%s.EIG" % filename, "r") as f:
        content=f.readlines()
        fermi_level=float(content[0].split()[0])
        n_bands= int(content[1].split()[0])
        n_spin = int(content[1].split()[1])
        n_kpoints = int(content[1].split()[2])
        eigenvalues=np.zeros((n_kpoints,n_spin,n_bands))
    
        energy=[]
        counter=1.0
        line_number=2
        for k in range(0,n_kpoints):
            energy+=content[line_number].split()[1:]
            for n in range(1, int(np.ceil(n_bands*n_spin/10.0))):
                line_number+=1
                energy+=content[line_number].split()
            for s in range (n_spin):
                for i in range(0, n_bands):
                    eigenvalues[k,s,i] = energy[i+n_bands*s]
            line_number+=1
            energy=[]
    
        return(fermi_level,eigenvalues)    


def siesta_read_kpoints(filename):
    """
    This routine reads the weights of K points and returns a 2 d array with dimensions:
    1) index of the k point 
    2) x,y,z coordinate and weight
    
    written by Mikhail Askerka, Yale University, 03/17/2015
    """

    with open("%s.KP" % filename, "r") as f:
        content=f.readlines()
        n_kpoints=int(content[0].split()[0])
        kpoints_weights=np.zeros((n_kpoints,4))
        line_number=1
        while line_number <=len(content)-1:
            kpoints_weights[line_number-1]= content[line_number].split()[1:] 
            line_number+=1
    #from 1/bohr to 1/angstrom
        kpoints_weights[:,:3] /= bohr
    return kpoints_weights


def siesta_read_struct_out(filename):
    """
    reads the .STRUCT_OUT file for the unit cell and positions
    is needed by siesta_read_HSX
    cell is a 3x3 matrix in Angstrom units
    
    written by Reinhard J. Maurer and Mikhail Askerka Yale University, 03/17/2015
    """
    with open("%s.STRUCT_OUT" % filename, "r") as f:
        content=f.readlines()
        cell = np.zeros([3,3],dtype=np.float)
        atomic_positions = np.zeros([len(content)-4,3],dtype=np.float)
        for i in range(3):
            bla = content[i].split()
            for j in range(3):
                cell[i,j] = np.float(bla[j])
        for i in range(len(content)-4):
            for j in range(3):
                atomic_positions[i,j] = np.float(content[i+4].split()[j+2])
	for i in range(len(atomic_positions)):
            atomic_positions[i]=np.dot(atomic_positions[i],cell)
    return cell, atomic_positions

def siesta_read_coefficients(filename, debug=0):
    """
    This routine reads siesta eigenvectors from MyM.WFSX binary
    using the FortranFile magical binary reader
    and returns an array with the following dimensions:
    1) the k-point index
    2) the spin index
    3) the state index (nwflist)
    4) the coefficient for each basis set index (nuotot)
    
    written by Mikhail Askerka, Yale University, 03/17/2015
    """

    f= FortranFile("%s.WFSX" % filename)
    nk,gamma= f.read_ints()
    if debug: print "nk, gamma  = ", nk,gamma
    nspin = int(f.read_ints())
    if debug: print "nspin =",nspin
    nuotot = int(f.read_ints())
    if debug: print "nuotot =",nuotot
    bla = f.read_record([('orbital_pos','i4'),('b','20S'),('c','i4'),('d','i4'),('e','20S')])
    orbital_pos = np.array(bla['orbital_pos']-1,dtype=np.int)
    psi=np.zeros((nk,nspin,nuotot,nuotot))
    psi = psi +0j
    #separate condition for gamma point since the array is real
    if bool(gamma)==True:  
        for iispin in range (1,nspin+1):  #for each spin
            f.read_record('f')
            ispin = int(f.read_ints())
            nwflist =int(f.read_ints())
            for iw in range(1,nwflist+1):  
                # for each state (nwflist = total number of states)
                indwf=f.read_ints()
                # we first read the energy of the state 
                energy=f.read_reals('d')
                # and all the orbital coefficients
                read_psi = f.read_reals('f')   
            read_psi=np.reshape(read_psi, (nuotot,1)) 
            psi[0,iispin-1,iw-1,:]=read_psi[:,0]
    else:
        for iik in range (1,nk+1):  #for each k-point
            for iispin in range (1,nspin+1):  #for each spin
                f.read_record('f')
                ispin = int(f.read_ints())
                nwflist =int(f.read_ints())
                for iw in range(1,nwflist+1):  
                    # for each state (nwflist = total number of states)
                    indwf=f.read_ints()
                    # we first read the energy of the state 
                    energy=f.read_reals('d')
                    # and all the orbital coefficients (real value, followed by the imaginary value
                    read_psi = f.read_reals('f')
                    if debug: print "nutot", nuotot, 'len(read_psi)', len(read_psi)
                    read_psi=np.reshape(read_psi, (nuotot,2))  # reshape it 
                    # and make a row of complex numbers
                    psi[iik-1,iispin-1,iw-1,:]=read_psi[:,0]+1j*read_psi[:,1]  
    return psi, orbital_pos

def siesta_read_HSX(kpts_array, filename, debug=0):
    """
    reads Hamiltonian, and Overlap from Siesta .HSX file

    Input: kpoints_array, filename (everything before .HSX)
    kpoints array must have dimensions [nkpts, 4] where the 
    first three elements are the fractional kpoint coordinates
    Output: H[kpoints,spin,norb, norb], S[kpoints, norb, norb]
    units are in Rydberg and Bohr!!
    
    written by Reinhard J. Maurer, Yale University, 03/17/2015
    """

    #two different cases
    #gamma=0: just gamma point (TODO not yet implemented) gamma!=0: k sampling
    f = FortranFile(filename+'.HSX')
    nkpts = len(kpts_array)
    #transform back to 1/bohr
    kpts = kpts_array.copy()
    kpts[:,:3] *= bohr
    if debug: print 'Reading Hamiltonian and Overlap from file : '+filename+'.HSX'
    
    #1st line: 
    # No. of orbs in cell, No. of orbs in supercell, nspin, nnzs
    no_u, no_s, nspin, nnzs = f.read_ints()
    if debug: print no_u, no_s, nspin, nnzs
    # 0 if klist, -1 if only gamma pnt is used
    gamma = bool(f.read_ints()[0])
    if debug: print "gamma from HSX", gamma
    #index list, for every orb in supercell finds 
    #correspond orb in small cell
    if gamma==False:
        indxuo = f.read_ints() 
        if debug: 
            print 'indxuo'
            print indxuo.shape
            print indxuo
    #numh contains the number of orbitals connected 
    #with orbital i
    numh = f.read_ints()
    if debug:
        print 'numh'
        print numh.shape
        print numh
    #generating listhptr
    listhptr = np.zeros(no_u,dtype=np.float)
    listhptr[0] = 0
    for oi in xrange(1,no_u):
        listhptr[oi] = listhptr[oi-1] + numh[oi-1]
    if debug:
        print 'listhptr'
        print listhptr
    #reading listh
    listh = np.zeros(nnzs,dtype=np.float)
    for oi in xrange(no_u):
        bla = f.read_ints()
        for im in xrange(numh[oi]):
            listh[listhptr[oi]+im] = bla[im] 
    if debug: 
        print 'listh'
        print listh.shape
        print listh

    h = np.zeros([nnzs,nspin],dtype=np.float)
    for si in range(nspin):
        for oi in range(no_u):
            bla = f.read_reals('f')
            for mi in xrange(numh[oi]):
                h[listhptr[oi]+mi,si] = bla[mi]
    s = np.zeros([nnzs],dtype=np.float)
    for oi in range(no_u):
        bla = f.read_reals('f')
        for mi in xrange(numh[oi]):
            s[listhptr[oi]+mi] = bla[mi]
    if debug:
        print h.shape
        print s.shape
   
    #reading No. of elecs, fermi smearing
    nelec, fermi_smear = f.read_reals('d')
    if debug: print nelec, fermi_smear

    if gamma==False:
        xij = np.zeros([nnzs,3],dtype=np.float)
        #now we read the xij array, atomic position distances in supercell
        for oi in range(no_u):
            bla = f.read_reals('f')
            for mi in xrange(numh[oi]):
                for xyz in range(3):
                    xij[listhptr[oi]+mi,xyz] = bla[(3*mi)+xyz]
    if debug: print xij.shape

    ##THAT was the READ_IN part
    ## now we construct the hamiltonian and overlap
    if gamma==False:
        H = np.zeros([nkpts,nspin,no_u,no_u],dtype=np.complex)
        S = np.zeros([nkpts,no_u,no_u],dtype=np.complex)
    else:
        H = np.zeros([1,nspin,no_u,no_u],dtype=np.float)
        S = np.zeros([1,no_u,no_u],dtype=np.float)

    if gamma==False:
        kpts = np.array(kpts,dtype=np.float)
        xij = np.array(xij,dtype=np.float)
        h = np.array(h,dtype=np.float)
        s = np.array(s,dtype=np.float)
        listh = np.array(listh,dtype=np.int)
        listhptr = np.array(listhptr,dtype=np.int)
        numh = np.array(numh,dtype=np.int)
        indxuo = np.array(indxuo,dtype=np.int)
        H_r, H_i, S_r, S_i = siesta_mod.siesta_calc_HSX(nspin, kpts, 
                no_u, numh, listhptr, listh, indxuo, xij, h, s)
        H[:,:,:,:] = H_r[:,:,:,:] +1.j*H_i[:,:,:,:]
        S[:,:,:] = S_r[:,:,:] +1.j*S_i[:,:,:]
    # only gamma point
    else: # TODO implement for gamma point only
        H[0,:,:,:] = h.reshape(1,nspin,no_u,no_u)
        S[0,:,:] = s.reshape(1,no_u,no_u)
    
    #this value of Rydberg is directly from siesta/Src/units.f90
    # Siesta Version 3.2 04/29/2015
    H *= complex(ryd)
    return H, S


def siesta_read_total_energy(filename):
    """
    Reads the total energy from the siesta output file 
    The way the utility is written now the file will be called 
    input.out
    """

    with open("%s.out" % filename, "r") as f:
        for line in f.readlines():
            if "siesta:         Total"  in line:
                total_energy = line.split()[3]

    return total_energy


"""
To call the routine from this file just type 
from siesta py import function
function(a,b)
"""
