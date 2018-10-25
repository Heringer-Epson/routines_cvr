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

    def __init__(self, case, dirpath, rest, step, ramp):
        
        self.dirpath = dirpath
        self.rest = rest
        self.step = step
        self.ramp = ramp
        
        #Default quantities are defined below.
        self.smoothing_window = 21

        #Default parameter space to explore:
        self.A = np.arange(0., 3., 0.05)
        self.tau = np.logspace(-2., 0., 20)
        self.B = np.arange(0., 150., 10.)
        
        self.parspace = np.asarray(
          [(_A,_tau,_B) for _A in self.A for _tau in self.tau for _B in self.B])

        self.parspace_reduced = np.asarray(
          [(_A,_tau) for _A in self.A for _tau in self.tau])

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
            self.show_fig = True
            self.save_fig = False
                
        create_run_dir('./../OUTPUT_FILES/RUNS/' + self.subdir)
            
