"""
Parsing routines for Siesta

Copyright Mikhail Askerka and Reinhard Maurer, Yale University, 03/17/2015
"""

import numpy as np
import struct
from scipy.io import FortranFile
import siesta_mod

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
        fermi_level=content[0].split()[0]
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
    return kpoints_weights


def siesta_read_struct_out(filename):
    """
    reads the .STRUCT_OUT file for the unit cell and positions
    is needed by siesta_read_HSX
    cell is a 3x3 matrix in Angstrom units
    
    written by Reinhard J. Maurer, Yale University, 03/17/2015
    """
    with open("%s.STRUCT_OUT" % filename, "r") as f:
        content=f.readlines()
        cell = np.zeros([3,3],dtype=np.float)
        for i in range(3):
            bla = content[i].split()
            for j in range(3):
                cell[i,j] = np.float(bla[j])
    
    return cell 
def siesta_read_coefficients(filename):
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
    print "nk, gamma  = ", nk,gamma
    nspin = int(f.read_ints())
    print "nspin =",nspin
    nuotot = int(f.read_ints())
    print "nuotot =",nuotot
    f.read_record('f')
    psi=np.zeros((nk,nspin,nuotot,nuotot))
    psi = psi +0j
    for iik in range (1,nk+1):  #for each k-point
        for iispin in range (1,nspin+1):  #for each spin
            f.read_record('f')
            ispin = int(f.read_ints())
            nwflist =int(f.read_ints())
            for iw in range(1,nwflist+1):  # for each state (nwflist = total number of states)
                indwf=f.read_ints()
                energy=f.read_reals('d')    # we first read the energy of the state 
                read_psi = f.read_reals('f')   # and all the orbital coefficients (real value, followed by the imaginary value
                read_psi=np.reshape(read_psi, (nwflist,2))  # reshape it 
                psi[iik-1,iispin-1,iw-1,:]=read_psi[:,0]+1j*read_psi[:,1]  # and make a row of complex numbers
    return psi

def siesta_read_HSX(kpts_array, cell, filename, debug=0):
    """
    reads Hamiltonian, and Overlap from Siesta .HSX file

    Input: kpoints_array, filename (everything before .HSX)
    kpoints array must have dimensions [nkpts, 4] where the 
    first three elements are the fractional kpoint coordinates
    Output: H[norb, norb, nkpts, nspin], S[norb, norb]
    units are in Rydberg and Bohr!!
    
    written by Reinhard J. Maurer, Yale University, 03/17/2015
    """

    #two different cases
    #just gamma point or k sampling
    f = FortranFile(filename+'.HSX')
    nkpts = len(kpts_array)

    print 'Reading Hamiltonian and Overlap from file : '+filename+'.HSX'
    
    #Angstrom to Bohr
    cell = cell*1.889726124565
    rcell = 2. * np.pi * np.linalg.inv(cell.T)
    ##transform kpts_array to absolute kpts
    for i, k in enumerate(kpts_array):
        kpts_array[i,:3] = np.dot(kpts_array[i,:3],rcell)

    #1st line: 
    # No. of orbs in cell, No. of orbs in supercell, nspin, nnzs
    no_u, no_s, nspin, nnzs = f.read_ints()
    if debug: print no_u, no_s, nspin, nnzs
    # 0 if klist, 1 if only gamma pnt is used
    gamma = f.read_ints()[0]
    if debug: print gamma
    #index list, for every orb in supercell finds 
    #correspond orb in small cell
    if gamma==0:
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

    if gamma==0:
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
    if gamma==0:
        H = np.zeros([nkpts,nspin,no_u,no_u],dtype=np.complex)
        S = np.zeros([nkpts,no_u,no_u],dtype=np.complex)
    else:
        H = np.zeros([1,nspin,no_u,no_u],dtype=np.float)
        S = np.zeros([1,no_u,no_u],dtype=np.float)

    # only gamma point
    if gamma==0:
        kpts_array = np.array(kpts_array,dtype=np.float)
        xij = np.array(xij,dtype=np.float)
        h = np.array(h,dtype=np.float)
        s = np.array(s,dtype=np.float)
        listh = np.array(listh,dtype=np.int)
        listhptr = np.array(listhptr,dtype=np.int)
        numh = np.array(numh,dtype=np.int)
        indxuo = np.array(indxuo,dtype=np.int)
        H_r, H_i, S_r, S_i = siesta_mod.siesta_calc_HSX(nspin, kpts_array, 
                no_u, numh, listhptr, listh, indxuo, xij, h, s)
        H[:,:,:,:] = H_r[:,:,:,:] +1j*H_i[:,:,:,:]
        S[:,:,:] = S_r[:,:,:] +1j*S_i[:,:,:]
    else:
        H[0,:,:,:] = h.reshape(nspin,no_u,no_u)
        S[0,:,:] = s.reshape(1,no_u,no_u)
    
    H *= complex(13.60569253)
    return H, S

"""
To call the routine from this file just type 
from siesta py import function
function(a,b)
"""
