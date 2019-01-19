import sys, os, time
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import cPickle
import numpy as np
from scipy.optimize import curve_fit
from multiprocessing import Pool
from functools import partial

def f(x, A, B):
    return A * x + B
def residual(y0,y1):
    return np.sum((y0 - y1)**2.)
def pars2line(idx, tau, A, B):
    return ('\n' + str(idx) + ',' + str(format(tau, '.6f'))
            + ',' + str(format(A, '.6f')) + ',' + str(format(B, '.6f')))

lbnd, ubnd = np.array([-10.,20.]), np.array([10.,200.]) #[Amin,Bmin] [Amax,Bmax]
def fit_line(taus, S, I, t, idx_space):
    min_res, best_tau, best_i = +np.inf, np.nan, np.nan
    bold_voxel = S['signal_ns'][idx_space,:]
    signal_unc = S['signal_noise'][idx_space,:]
    for i, tau in enumerate(taus):
        pCO2_conv = I(t, tau)
        try:
            est, pcov = curve_fit(
              f, pCO2_conv, bold_voxel, p0=[0.5,70.],
              sigma=signal_unc, bounds=(lbnd, ubnd))
            res = residual(bold_voxel, f(pCO2_conv, est[0], est[1]))
            if res <= min_res:
                best_tau, best_i = tau, i
                best_A, best_B  = est[0], est[1]
                min_res = res
        except:
            est = (np.nan, np.nan, np.nan)
            res = np.inf
    return pars2line(idx_space,best_tau,best_A,best_B)

class Estimate_Bestpars(object):
    """
    Code Description
    ----------    
    Estimate the most likely parameters (A,B,tau), which will then be used
    in a MCMC routine.
    
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
        fpath_out = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'estimated_A_tau_B.csv'
        with open(fpath_S, 'r') as fS, open(fpath_I, 'r') as fI,\
          open(fpath_out, 'w') as out:
            out.write('voxel,tau_est,A_est,B_est')
            S, I = cPickle.load(fS), cPickle.load(fI)
            t = S['time']
            N_voxels = S['signal_ns'].shape[0]
            start_time = time.time()
            
            output = []
            fit_given_idx = partial(fit_line, _run.tau, S, I, t)
            pool = Pool(5)
            output += pool.map(fit_given_idx,range(N_voxels))
            #output += pool.map(fit_given_idx,range(10))
            pool.close()
            pool.join()
            for line in output:
                out.write(line) 
            delta_time = time.time() - start_time

        print '    Run took ', format(delta_time, '.1f'), 's'
        print '    Approx ', format(delta_time / float(N_voxels), '.3f'), ' s/voxel'
        print '    Done.\n'
