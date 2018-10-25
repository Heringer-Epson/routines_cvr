#!/usr/bin/env python
import os
import numpy as np
import cPickle
from scipy.integrate import simps
from scipy.interpolate import interp1d
from multiprocessing import Pool
from functools import partial

def exp_func(t, tau):
    return np.exp(-t / tau)

def compute_convolution(func_a,integration_range,tau,t):
    cond = (integration_range <= t)        
    integration_points = integration_range[cond]
    a = func_a(integration_points)
    b = exp_func(t - integration_points, tau)
    _expected_signal = simps(np.multiply(a,b), integration_points)
    return _expected_signal

def Compute_Model(_pCO2, _time, _ts, c_old, _A, _tau, _B):
    conv_pCO2 = []
    _time_fine = np.arange(_time[0], _time[-1], _ts)        

    pCO2_func = interp1d(_time,_pCO2,bounds_error=False,
                         fill_value=(_pCO2[0],_pCO2[1]))

    for t_inp in _time:
        conv_pCO2.append(compute_convolution(pCO2_func, _time_fine, _tau, t_inp))

    #Ensure that the convolved function has the same area as the
    #input pCO2 function. This is important in derivinf the ssCVR.
    c_new = simps(conv_pCO2,_time)
    conv_pCO2 = np.array(conv_pCO2) * (_A * c_old / c_new) + _B
    
    return _A, _tau, _b, conv_pCO2

class Save_Models(object):
    def __init__(self, _run):
        M = {}
        M['parspace'] = _run.parspace
        output = []


        c_old = simps(_run.pCO2,_run.time) 

        N_sets = len(_run.A)
        N_subsets = len(_run.tau)
        for i, A in enumerate(_run.A):
            print 'Calculating set ' + str(i + 1) + '/' + str(N_sets)
            for j, tau in enumerate(_run.tau):
                print '    Calculating set ' + str(j + 1) + '/' + str(N_subsets)
                            
                L_of_tau = partial(Compute_Model,_run.pCO2,_run.time,0.1,c_old,A,tau)

                pool = Pool(5)
                output += pool.map(L_of_tau,_run.B)            
                pool.close()
                pool.join()       
        
        #It's been tested that the output order matches the input one.
        for p, conv in enumerate(zip(*output)[3]):
            label = 'conv_' + str(p)
            M[label] = conv
        fpath = os.path.join(
          './../OUTPUT_FILES/RUNS/' + _run.subdir, 'models.pkl')
        with open(fpath, 'w') as output:
            cPickle.dump(M, output, cPickle.HIGHEST_PROTOCOL)
        print 'Done'

'''

#!/usr/bin/env python
import os
import numpy as np
import cPickle
from scipy.integrate import simps
from scipy.interpolate import interp1d
from multiprocessing import Pool
from functools import partial

def exp_func(t, tau):
    return np.exp(-t / tau)

def compute_convolution(func_a,integration_range,tau,t):
    cond = (integration_range <= t)        
    integration_points = integration_range[cond]
    a = func_a(integration_points)
    b = exp_func(t - integration_points, tau)
    _expected_signal = simps(np.multiply(a,b), integration_points)
    return _expected_signal

def Compute_Model(_pCO2, _time, _ts, c_old, _A, _tau):
    conv_pCO2 = []
    _time_fine = np.arange(_time[0], _time[-1], _ts)        

    pCO2_func = interp1d(_time,_pCO2,bounds_error=False,
                         fill_value=(_pCO2[0],_pCO2[1]))

    for t_inp in _time:
        conv_pCO2.append(compute_convolution(pCO2_func, _time_fine, _tau, t_inp))

    #Ensure that the convolved function has the same area as the
    #input pCO2 function. This is important in derivinf the ssCVR.
    c_new = simps(conv_pCO2,_time)
    conv_pCO2 = np.array(conv_pCO2) * (_A * c_old / c_new)
    
    return _A, _tau, conv_pCO2

class Save_Models(object):
    def __init__(self, _run):
        M = {}
        M['parspace'] = _run.parspace
        output = []


        c_old = simps(_run.pCO2,_run.time) 

        N_sets = len(_run.A)
        for i, A in enumerate(_run.A):
            print 'Calculating set ' + str(i + 1) + '/' + str(N_sets)
                        
            L_of_tau = partial(Compute_Model,_run.pCO2,_run.time,0.1,c_old,A)

            pool = Pool(5)
            output += pool.map(L_of_tau,_run.tau)            
            pool.close()
            pool.join()       
        
        #It's been tested that the output order matches the input one.
        for p, conv in enumerate(zip(*output)[2]):
            label = 'conv_' + str(p)
            M[label] = conv
        fpath = os.path.join(
          './../OUTPUT_FILES/RUNS/' + _run.subdir, 'models.pkl')
        with open(fpath, 'w') as output:
            cPickle.dump(M, output, cPickle.HIGHEST_PROTOCOL)
        print 'Done'
'''
