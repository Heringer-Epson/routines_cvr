import time
import cPickle
import numpy as np
from scipy.optimize import curve_fit
from multiprocessing import Pool
from functools import partial

def f(x, A, B):
    return A * x + B
def residual(y0,y1):
    return np.sum((y0 - y1)**2.)
def pars2line(idx, tau, A, A_unc, B, B_unc):
    return ('\n' + str(idx) + ',' + str(format(tau, '.6f'))
            + ',' + str(format(A, '.6f')) + ',' + str(format(A_unc, '.6f')) 
            + ',' + str(format(B, '.6f')) + ',' + str(format(B_unc, '.6f')))

#lbnd, ubnd = np.array([-10.,20.]), np.array([10.,200.]) #[Amin,Bmin] [Amax,Bmax]
def fit_line(taus, S, I, t, bound_low, bound_up, idx_space):
    min_res, best_tau, best_i = +np.inf, np.nan, np.nan,
    best_A, best_B, perr = np.nan, np.nan, (np.nan, np.nan)
    bold_voxel = S['signal_ns'][idx_space,:]
    signal_unc = S['signal_noise'][idx_space,:]
    for i, tau in enumerate(taus):
        pCO2_conv = I(t, tau)
        try:
            est, pcov = curve_fit(
              f, pCO2_conv, bold_voxel, p0=[0.2,90.],
              sigma=signal_unc, bounds=(bound_low, bound_up))
            perr = np.sqrt(np.diag(pcov))
            res = residual(bold_voxel, f(pCO2_conv, est[0], est[1]))
            if res <= min_res:
                best_tau, best_i = tau, i
                best_A, best_B  = est[0], est[1]
                min_res = res
        except:
            est = (np.nan, np.nan, np.nan)
            res = np.inf
    return pars2line(idx_space,best_tau,best_A,perr[0],best_B,perr[1])

class Estimate_Bestpars(object):
    """
    Code Description
    ----------    
    Estimate the most likely parameters (A,B,tau), which will then be used
    in a MCMC routine. Because tau is non-linear, it is estiamted in discrete
    steps, while A and B are estimated using curve_fit, under a prior boundary.
    The parameter 'C' (linear correction following the time i.e. C*t) is
    not done because writing the signal function with both pCO2_conv and t
    would be too intricate. Instead, C will always be guessed as zero when
    running the MCMC routine.
    
    Parameters
    ----------
    _run : object
      instance of input_params containing the data in each case of study.

    Outputs:
    --------
    ./../OUTPUT_FILES/RUNS/' + self._run.subdir + 'estimated_A_tau_B.csv'     
    """
    #@profile
    def __init__(self, _run):
        print 'Estimating parameters...'       
        fpath_S = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'PICKLES/smooth.pkl'
        fpath_I = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'PICKLES/tau_interp.pkl'
        fpath_out = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'estimated_pars.csv'
        with open(fpath_S, 'r') as fS, open(fpath_I, 'r') as fI,\
          open(fpath_out, 'w') as out:
            out.write('voxel,tau,A,A_unc,B,B_unc')
            S, I = cPickle.load(fS), cPickle.load(fI)
            t = S['time']
            N_voxels = S['signal_ns'].shape[0]
            start_time = time.time()

            bound_low = np.array([_run.A_lim[0],_run.B_lim[0]])
            bound_up = np.array([_run.A_lim[1],_run.B_lim[1]])
                    
            output = []
            fit_given_idx = partial(
              fit_line, _run.tau, S, I, t, bound_low, bound_up)
            pool = Pool(5)
            output += pool.map(fit_given_idx,range(N_voxels))
            #output += pool.map(fit_given_idx,range(100))
            pool.close()
            pool.join()
            for line in output:
                out.write(line) 
            
        print '    Run took ', format(time.time() - start_time, '.1f'), 's'
        print '    Done.\n'
