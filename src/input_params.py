#!/usr/bin/env python
import os
import numpy as np

def create_run_dir(output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)  
        os.mkdir(output_dir + 'FIGURES/')  

class Input_Parameters(object):
    """Set input parameters for CVR analysis.
    """

    def __init__(self, case, dirpath, rest, step, ramp, show_fig, save_fig):
        
        self.dirpath = dirpath
        self.rest = rest
        self.step = step
        self.ramp = ramp
        self.show_fig = show_fig
        self.save_fig = save_fig
        
        #Default quantities are defined below.
        self.smoothing_window = 21

        #Default parameter space to explore:
        self.A = np.arange(-1., 1., 0.05)
        #self.tau = np.arange(.0, 8.4, 0.05)
        self.tau = np.logspace(-2., 10., 10)
        self.B = np.arange(0., 150., 10.)
        self.parspace = np.asarray(
          [(_A,_tau,_B) for _A in self.A for _tau in self.tau for _B in self.B])
        
        #Number of baselines to marginalize over.
        
        #self.B_trials = 10
        #self.B = np.zeros(self.B_trials)
        
        #Initialize quantities to be compute.
        self.signal = None
        self.pCO2 = None
        self.time = None
        self.M = None

        self._signal_n, self.signal_sn, self.signal_ns_unc = None, None, None
        self._pCO2_s, self.pCO2_noise = None, None      
        
        #Used to compute the default tomography analysis.
        if case == 'default':
            self.subdir = 'test1/'    
            self.time_step = 2.4 #in units of second.
            
        create_run_dir('./../OUTPUT_FILES/RUNS/' + self.subdir)
            
