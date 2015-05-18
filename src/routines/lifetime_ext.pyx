cimport cython
import numpy as np
cimport numpy as np
REAL_TYPE = np.float
INT_TYPE = np.int
COMPLEX_TYPE = np.complex
ctypedef np.int_t INT_TYPE_t
ctypedef np.float_t REAL_TYPE_t
ctypedef np.complex_t COMPLEX_TYPE_t

from libc.math cimport sqrt, exp, cos, sin
from vibcooling.constants import hbar_1, pi, sqrt_pi, conversion_factor

def calculate_nonadiabatic_coupling_HS_kpts(
        int state_i, int state_f,
        np.ndarray[COMPLEX_TYPE_t, ndim=3] psi, 
        np.ndarray[REAL_TYPE_t, ndim=1] Nvec1, 
        np.ndarray[REAL_TYPE_t, ndim=1] Nvec2, 
        np.ndarray[COMPLEX_TYPE_t, ndim=3] H_plus, 
        np.ndarray[COMPLEX_TYPE_t, ndim=3] H_minus, 
        np.ndarray[COMPLEX_TYPE_t, ndim=3] S_plus, 
        np.ndarray[COMPLEX_TYPE_t, ndim=3] S_minus, 
        np.ndarray[REAL_TYPE_t, ndim=2] atomic_positions, 
        np.ndarray[INT_TYPE_t,  ndim=1] basis_pos, 
        np.ndarray[REAL_TYPE_t, ndim=2] kpts, 
        np.ndarray[REAL_TYPE_t, ndim=1] kweights,
        float dq, 
        float fermi_level, 
        float e,
        ):
    """
    Cython extension to calculate a general nonadiabatic coupling element in the 
    H-eS formalism of Head-Gordon where states can be at differen k points.
    """
  
    cdef Py_ssize_t a, b, l, k 
    cdef INT_TYPE_t n_basis, nk
    cdef REAL_TYPE_t product_real,product_r, product_i, tmp_r, tmp_i, kw
    cdef REAL_TYPE_t kx, H_q_r,H_q_i,S_q_r,S_q_i
    cdef np.ndarray[REAL_TYPE_t,ndim=1] kvec = np.empty(3,dtype=REAL_TYPE)
    cdef np.ndarray[REAL_TYPE_t,ndim=1] RNvec1 = np.empty(3,dtype=REAL_TYPE)
    cdef np.ndarray[REAL_TYPE_t,ndim=1] RNvec2 = np.empty(3,dtype=REAL_TYPE)

    nbasis = len(basis_pos) 
    product_r = 0.0
    product_i = 0.0
    nk = len(kpts)
    for k in range(nk):
        kvec[:] = kpts[k]
        kw = kweights[k]
        for a in range(nbasis):
            RNvec1[:] = atomic_positions[basis_pos[a]]+Nvec1[:]
            for b in range(nbasis):
                RNvec2[:] = atomic_positions[basis_pos[b]]+Nvec2[:]
                tmp_r = 0.0
                tmp_i = 0.0
                #Fourier transformation of G matrix element
                H_q_r=(H_plus[k,a,b].real-H_minus[k,a,b].real)/2/dq
                H_q_i=(H_plus[k,a,b].imag-H_minus[k,a,b].imag)/2/dq
                S_q_r=(S_plus[k,a,b].real-S_minus[k,a,b].real)/2/dq
                S_q_i=(S_plus[k,a,b].imag-S_minus[k,a,b].imag)/2/dq
                kx = kvec[0]*(RNvec1[0]-RNvec2[0])+ \
                     kvec[1]*(RNvec1[1]-RNvec2[1])+ \
                     kvec[2]*(RNvec1[2]-RNvec2[2])
                tmp_r += cos(kx)*(H_q_r-fermi_level*S_q_r)
                tmp_i += sin(kx)*(H_q_i-fermi_level*S_q_i)

                #wavefunc coeffs
                tmp_r *= psi[k,state_i,a].real*psi[k,state_f,b].real
                tmp_i *= -psi[k,state_i,a].imag*psi[k,state_f,b].imag
                product_r += tmp_r
                product_i += tmp_i
        product_r *= kw
        product_i *= kw
    product_real = product_r*product_r + product_i*product_i 
    product_real /= (e*e)
    return product_real





