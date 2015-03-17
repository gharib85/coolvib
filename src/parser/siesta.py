"""
Parsing routine for Siesta
"""

import numpy as np
import struct
from scipy.io import FortranFile
import siesta_mod

def siesta_read_eigenvalues(filename):
    """
    This routine reads the number of K points, the eigenvalues and the \
    and the number of states

    Input: the name of the eigenvalue file
    Output: Fermi_level, Eigenvalues i

    The Eigenvalues is a 3D numpy array:
    1) K point
    2) Spin
    3) State
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


def siesta_read_coefficients(filename):
    """
    This routine reads siesta eigenvectors from MyM.WFSX binary
    using the FortranFile magical binary reader
    and returns an array with the following dimensions:
    1) the k-point index
    2) the spin index
    3) the state index (nwflist)
    4) the coefficient for each basis set index (nuotot)
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

def siesta_read_HSX(kpts_array, filename, debug=0):
    """
    reads Hamiltonian, and Overlap from Siesta .HSX file

    Input: kpoints_array, filename (everything before .HSX)
    kpoints array must have dimensions [nkpts, 4] where the 
    first three elements are the fractional kpoint coordinates
    Output: H[norb, norb, nkpts, nspin], S[norb, norb]
    units are in Rydberg!!
    """

    #two different cases
    #just gamma point or k sampling
    f = FortranFile(filename+'.HSX')
    nkpts = len(kpts_array)
    
    print 'Reading Hamiltonian and Overlap from file : '+filename+'.HSX'

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
    #listhptr = np.zeros(nu_o,dtype=np.float)
    #listhptr[0] = 0
    #for oi in xrange(1,nu_o):
        #listhptr[oi] = listhptr[oi-1] + numh[oi-1]
    #listh = np.zeros(nnzs,dtype=np.float)
    #for oi in xrange(nu_o):
        #listh[listhptr[oi]:(listhptr[oi]+numh[oi])] = f.read_ints()
    listh  = []
    for i in xrange(no_u):
        listh += list(f.read_ints())
    listh = np.array(listh)
    if debug: 
        print 'listh'
        print listh.shape
        print listh

    h = np.zeros([nnzs,nspin],dtype=np.float)
    for si in range(nspin):
        tmp = []
        for i in range(no_u):
            bla = list(f.read_reals('f'))
            tmp += bla 
        h[:,si] = np.array(tmp)
    s = []
    for i in range(no_u):
        bla = list(f.read_reals('f'))
        s += bla
    s = np.array(s)
    if debug:
        print h.shape
        print s.shape
   
    #reading No. of elecs, fermi smearing
    nelec, fermi_smear = f.read_reals('d')
    if debug: print nelec, fermi_smear

    if gamma==0:
        #now we read the xij array, atomic position distances in supercell
        xij = []
        for i in range(no_u):
            xij += list(f.read_reals('f'))
    xij = np.array(xij)
    if debug: print xij[0].shape

    ##THAT was the READ_IN part
    ## now we construct the hamiltonian and overlap
    if gamma==0:
        H = np.zeros([no_u,no_u,nkpts,nspin],dtype=np.float)
        S = np.zeros([no_u,no_u,nkpts],dtype=np.float)
    else:
        H = np.zeros([no_u,no_u,1,nspin],dtype=np.float)
        S = np.zeros([no_u,no_u,1],dtype=np.float)
    
    # only gamma point
    if gamma==0:
        kpts_array = np.array(kpts_array,dtype=np.float)
        xij = np.array(xij,dtype=np.float)
        h = np.array(h,dtype=np.float)
        s = np.array(s,dtype=np.float)
        listh = np.array(listh,dtype=np.int)
        numh = np.array(numh,dtype=np.int)
        indxuo = np.array(indxuo,dtype=np.int)

        H_r, H_i, S_r, S_i = siesta_mod.siesta_calc_HSX(nspin, kpts_array, 
                no_u, numh, listh, indxuo, xij, h, s)
        H = H_r +1j*H_i
        S = S_r +1j*S_i

    else:
        H[:,:,0,:] = h.reshape(no_u,no_u,nspin)
        S[:,:] = s.reshape(no_u,no_u)

    return H, S

"""
To call the routine from this file just type 
from siesta py import function
function(a,b)
"""
