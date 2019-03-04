#!/usr/bin/env python
import sys, os, time
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import numpy as np
import cPickle
from scipy.integrate import simps
from scipy.interpolate import interp1d, interp2d

def exp_func(t, tau):
    return np.exp(-t / tau)

def compute_convolution(func_a,integration_range,tau,t):
    cond = (integration_range <= t)        
    integration_points = integration_range[cond]
    a = func_a(integration_points)
    b = exp_func(t - integration_points, tau)
    _expected_signal = simps(np.multiply(a,b), integration_points)
    return _expected_signal

def Compute_Model(_t, _t_fine, _pCO2_func, _c_old, _tau):
    conv_pCO2 = []
    for t_inp in _t:
        conv_pCO2.append(compute_convolution(_pCO2_func, _t_fine, _tau, t_inp))
    conv_pCO2 = np.array(conv_pCO2)
    cond = ((_t > 0.))    
    c_new = simps(conv_pCO2[cond],_t[cond]) #Preserve area of unconvolved points.
    
    conv_pCO2 = np.array(conv_pCO2) * (_c_old / c_new)
    return conv_pCO2

class Tau_Interpolator(object):
    def __init__(self, _run):

        self._run = _run
        self.run_routine()

    def retrieve_data(self):
        fpath_D = './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'PICKLES/data.pkl'
        fD = open(fpath_D , 'r')
        self.D = cPickle.load(fD)
        fD.close()

    def construct_aux_pCO2(self):
        """This routine artificially extends the rest period of the measured
        pCO2. The reason for this is that when making a convolution, the
        convolved pCO2 at t~0s will have been impacted by the pCO2 level prior
        to that, which is not recorded. This issue is particular significant
        if tau is large. To amend this issue, the pCO2 level prior to 0s is
        assumed to be the median of the pCO2 during the rest phase.
        t_pre is assumed 300s, which can account for an e**3 folding time for
        the extreme case of tau=100s. Note that the normalization has to be
        over t>0 only!
        """
        t = self.D['time']
        pCO2 = self.D['pCO2']
        t_pre = -5000.
        rest_cond = ((t >= self._run.rest[0]) & (t <= self._run.rest[1]))
        rest_pCO2 = np.median(pCO2[rest_cond])

        t_neg = np.arange(t_pre,0.,self._run.time_step)
        pCO2_neg = np.zeros(len(t_neg)) + rest_pCO2
        
        self.t_aux = np.concatenate((t_neg,t),axis=0)
        self.pCO2_aux = np.concatenate((pCO2_neg,pCO2),axis=0)
                
    def interpolate_convolved_pCO2(self):

        print 'Creating matrix of convolved PCO2...'
        start_time = time.time()
        M = {}
        fpath_D = './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'PICKLES/data.pkl'
        fpath_out = './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'PICKLES/tau_interp.pkl'
        with open(fpath_out, 'w') as fO:
            t_step = 0.1
            t_fine = np.arange(self.t_aux[0], self.t_aux[-1], t_step)  
            pCO2_func = interp1d(
              self.t_aux,self.pCO2_aux,bounds_error=False,
              fill_value=(self.pCO2_aux[0],self.pCO2_aux[-1]))
            c_old = simps(self.D['pCO2'],self.D['time'])
            xx, yy = np.meshgrid(self.t_aux, self._run.tau)
            conv_pCO2_matrix = np.zeros(shape=(len(self._run.tau),len(self.t_aux)))
            for i, tau in enumerate(self._run.tau):
                conv_pCO2_matrix[i] = Compute_Model(self.t_aux, t_fine, pCO2_func, c_old, tau)
            f = interp2d(self.t_aux, self._run.tau, conv_pCO2_matrix, kind='cubic')
            cPickle.dump(f, fO)          
                  
        print '    Run took ', format(time.time() - start_time, '.1f'), 's'
        print '    Done.\n'

    def run_routine(self):
        self.retrieve_data()
        self.construct_aux_pCO2()
        self.interpolate_convolved_pCO2()
