#cython: boundscheck=False, wraparound=False, nonecheck=False
from libc.math cimport exp, log, log10, pow

cpdef double lnlike(double A, double B, double C, double[:] ts,
                    double[:] x_conv, double[:] y, double[:] yerr):
    cdef int N = y.shape[0]
    cdef double lnL = 0.
    cdef double model = 0.
    for i in range(N):
        model = A * x_conv[i] + B + C * ts[i]
        lnL += -pow(y[i]-model,2) / pow(yerr[i],2)
    return lnL

cpdef double chi2(double A, double B, double C, double[:] ts,
                    double[:] x_conv, double[:] y, double[:] yerr):
    cdef int N = y.shape[0]
    cdef double chi2 = 0.
    cdef double model = 0.
    for i in range(N):
        model = A * x_conv[i] + B + C * ts[i]
        chi2 += pow(y[i]-model,2) / pow(yerr[i],2)
    return chi2

cpdef double lnlike_gmwm(double A_gm, double A_wm, double[:] x_conv_gm,
                         double[:] x_conv_wm, double B, double C, 
                         double[:] ts, double[:] gf, double[:] wf, double[:] y,
                         double[:] yerr):
    cdef int N = y.shape[0]
    cdef double lnL = 0.
    cdef double model = 0.
    for i in range(N):
        model = gf[i] * A_gm * x_conv_gm[i] + wf[i] * A_wm * x_conv_wm[i] + B + C * ts[i]
        lnL += -pow(y[i]-model,2) / pow(yerr[i],2)
    return lnL

cpdef double chi2_gmwm(double A_gm, double A_wm, double B, double C, 
                       double[:] ts, double[:] gf, double[:] wf, double[:] x_conv_gm, double[:] x_conv_wm, 
                       double[:] y, double[:] yerr):
    cdef int N = y.shape[0]
    cdef double chi2 = 0.
    cdef double model = 0.
    for i in range(N):
        model = gf[i] * A_gm * x_conv_gm[i] + wf[i] * A_wm * x_conv_wm[i] + B + C * ts[i]
        chi2 += pow(y[i]-model,2) / pow(yerr[i],2)
    return chi2
