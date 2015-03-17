import cython
import numpy as np
cimport numpy as np
REAL_TYPE = np.float
INT_TYPE = np.int
COMPLEX_TYPE = np.complex
ctypedef np.int_t INT_TYPE_t
ctypedef np.float_t REAL_TYPE_t
ctypedef np.complex_t COMPLEX_TYPE_t

def siesta_calc_HSX(int nspin, np.ndarray[REAL_TYPE_t,ndim=2] kpts_array, 
                    int no_u, np.ndarray[INT_TYPE_t,ndim=1] numh, 
                    np.ndarray[INT_TYPE_t,ndim=1] listh, 
                    np.ndarray[INT_TYPE_t,ndim=1] indxuo,
                    np.ndarray[REAL_TYPE_t,ndim=1] xij,
                    np.ndarray[REAL_TYPE_t,ndim=2] h,
                    np.ndarray[REAL_TYPE_t,ndim=1] s,
            ):
    """
    calculates actual H and S matrices
    """

    #INIT H and S
    cdef int nkpts
    nkpts = len(kpts_array)
    cdef np.ndarray[REAL_TYPE_t,ndim=4] H_r = np.empty([no_u,no_u,nkpts,nspin],dtype=REAL_TYPE)
    cdef np.ndarray[REAL_TYPE_t,ndim=3] S_r = np.empty([no_u,no_u,nkpts],dtype=REAL_TYPE)
    cdef np.ndarray[REAL_TYPE_t,ndim=4] H_i = np.empty([no_u,no_u,nkpts,nspin],dtype=REAL_TYPE)
    cdef np.ndarray[REAL_TYPE_t,ndim=3] S_i = np.empty([no_u,no_u,nkpts],dtype=REAL_TYPE)

    cdef unsigned int si, ki, iuo, j, jo, juo
    cdef np.ndarray[REAL_TYPE_t,ndim=1] k=np.empty(4,dtype=REAL_TYPE), kvec=np.empty(3,dtype=REAL_TYPE)
    cdef REAL_TYPE_t phasef_r, phasef_i

    # only gamma point
    for si in range(nspin):
        for ki,k in enumerate(kpts_array):
            kvec = k[:3]
            for iuo in xrange(no_u):
                for j in xrange(numh[iuo]):
                    jo = listh[numh[iuo]+j]
                    #print 'jo :', jo
                    juo = indxuo[numh[iuo]+j] -1
                    #print 'juo :', juo
                    #print kvec
                    #print xij[iuo][j]
                    kx = np.dot(kvec,xij[(numh[iuo]+j)*3:(numh[iuo]+j)*3+3])
                    phasef_r = np.cos(kx)
                    phasef_i = np.sin(kx)
                    #print phasef
                    #print h[numh[iuo]+j,si]
                    #print iuo, juo, ki, si
                    H_r[iuo,juo,ki,si] += phasef_r*h[numh[iuo]+j,si]
                    H_i[iuo,juo,ki,si] += phasef_i*h[numh[iuo]+j,si]
                    if si==0:
                        S_r[iuo,juo,ki] += phasef_r*s[numh[iuo]+j]
                        S_i[iuo,juo,ki] += phasef_i*s[numh[iuo]+j]

    return H_r, H_i, S_r, S_i

