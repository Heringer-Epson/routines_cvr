#cython: boundscheck=False, wraparound=False, nonecheck=False
from libc.math cimport exp, log, log10, pow

cpdef double lnlike(double A, double B, double C, double D, double[:] ts,
                    double[:] x_conv, double[:] y, double[:] yerr):
    cdef int N = y.shape[0]
    cdef double lnL = 0.
    cdef double model = 0.
    for i in range(N):
        model = A * x_conv[i] + B + C * ts[i] + D * pow(ts[i],2)
        lnL += -pow(y[i]-model,2) / pow(yerr[i],2)
    return lnL

cpdef double chi2(double A, double B, double C, double D, double[:] ts,
                    double[:] x_conv, double[:] y, double[:] yerr):
    cdef int N = y.shape[0]
    cdef double chi2 = 0.
    cdef double model = 0.
    for i in range(N):
        model = A * x_conv[i] + B + C * ts[i] + D * pow(ts[i],2)
        chi2 += pow(y[i]-model,2) / pow(yerr[i],2)
    return chi2
