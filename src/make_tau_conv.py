#!/usr/bin/env python
import sys, os, time
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import numpy as np
import cPickle
from scipy.signal import savgol_filter
from scipy.integrate import simps
from scipy.interpolate import interp1d
from multiprocessing import Pool
from functools import partial
import data_handling 

def exp_func(t, tau):
    return np.exp(-t / tau)

def compute_convolution(func_a,integration_range,tau,t):
    cond = (integration_range <= t)        
    integration_points = integration_range[cond]
    a = func_a(integration_points)
    b = exp_func(t - integration_points, tau)
    _expected_signal = simps(np.multiply(a,b), integration_points)
    return _expected_signal

def Compute_Model(_pCO2, _time, _ts, c_old, _tau):
    conv_pCO2 = []
    _time_fine = np.arange(_time[0], _time[-1], _ts)        

    pCO2_func = interp1d(_time,_pCO2,bounds_error=False,
                         fill_value=(_pCO2[0],_pCO2[-1]))

    for t_inp in _time:
        conv_pCO2.append(compute_convolution(pCO2_func, _time_fine, _tau, t_inp))

    #Ensure that the convolved function has the same area as the
    #input pCO2 function. This is important in deriving the ssCVR.
    c_new = simps(conv_pCO2,_time)
    conv_pCO2 = np.array(conv_pCO2) * (c_old / c_new)
    return _tau, conv_pCO2

class Save_Models(object):
    def __init__(self, _run):
        print 'Calculating models...'
        M = {}
        output = []
        fpath_D = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'PICKLES/tau_conv.pkl'
        with open(fpath_D, 'r') as fD:
            D = cPickle.load(fD)
        
            c_old = simps(D['pCO2'],D['time']) 
            start_time = time.time()
            L_of_tau = partial(Compute_Model,D['pCO2'],D['time'],0.1,c_old)
            pool = Pool(5)
            output += pool.map(L_of_tau,_run.tau)            
            pool.close()
            pool.join()       
            delta_time = time.time() - start_time
            
            #For each tau, compute a model prediction depending also on A and B.
            #This avoids running repeatd convolution calculations.
            for (_tau,_conv_pCO2) in output:
                M[data_handling.pars2label(_A,_tau,_B)] = _A * _conv_pCO2 + _B
            data_handling.save_pickle(_run.subdir, 'models.pkl', M)
      
        print '    Run took ', format(delta_time, '.1f'), 's'
        print '    Approx ', format(delta_time / _run.N_cells * 1.e6, '.3f'), ' mus/cell'
        print '    Done.\n'
