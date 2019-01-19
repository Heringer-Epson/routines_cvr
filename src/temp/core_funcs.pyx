#cython: boundscheck=False, wraparound=False, nonecheck=False
from libc.math cimport exp, log, log10, pow

cdef extern from "vfastexp.h":
    double exp_approx "EXP" (double)

cpdef double compute_L(double[:] signal, double[:] model, double[:] unc):
    cdef double M_pi = 3.141592653
    cdef int N = signal.shape[0]
    cdef double ln_L = 0. 
    cdef int i
    for i in range(N):
        ln_L += (
          + log(1. / (2. * M_pi * pow(unc[i],2)))
          - pow(signal[i] - model[i],2) / (2. * pow(unc[i],2)))
    return ln_L

cpdef double marg(double[:] ln_L, double[:] theta):
    cdef int N = ln_L.shape[0]
    cdef double I = 0
    for i in range(N - 1): #Minus one because we need differences.
        #L_mean = (exp_approx(ln_L[i]) + exp_approx(ln_L[i + 1])) / 2.
        L_mean = (exp(ln_L[i]) + exp(ln_L[i + 1])) / 2.
        delta_theta = theta[i + 1] - theta[i]
        I += (L_mean * delta_theta)
    return log(I)
    
#=-=-=-=-=-=-=-=-=-=-=-=-=-=Functions used for MCMC-=-=-=-=-=-=-=-=-=-=-=-=-=-=


cpdef double lnlike(double A, double B, double[:] x_conv, double[:] y, double[:] yerr):
    cdef int N = y.shape[0]
    cdef double lnL = 0.
    cdef double model = 0.
    cdef double inv_sigma2 = 0.

    for i in range(N):
        model = A * x_conv[i] + B
        inv_sigma2 = 1.0 / (pow(yerr[i],2) + pow(model,2))
        lnL += pow(y[i]-model,2) * inv_sigma2 - log(inv_sigma2)
    return -0.5 * lnL



