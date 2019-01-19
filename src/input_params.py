#!/usr/bin/env python
import os
import numpy as np

def create_run_dir(output_dir):
    if not os.path.exists(output_dir):
        os.mkdir(output_dir)  
        os.mkdir(output_dir + 'FIGURES/')  
        os.mkdir(output_dir + 'PICKLES/')  

class Input_Parameters(object):
    """Set input parameters for CVR analysis.
    """

    def __init__(self, case, dirpath, rest, norest, step, ramp):
        print '\n\n\nSetting input parameters...'
        
        self.dirpath = dirpath
        self.rest = rest
        self.norest = norest
        self.step = step
        self.ramp = ramp
        
        #Default quantities are defined below.
        self.smoothing_window = 21

        
        self.A_step = 0.05
        self.A = np.arange(-1., 1., self.A_step)
        self.tau = np.logspace(0., 2., 100)
        self.B = np.arange(78., 82., .1)
        self.N_cells = float(len(self.A) * len(self.tau) * len(self.B))
        
        self.parspace = np.asarray(
          [(_A,_tau,_B) for _A in self.A for _tau in self.tau for _B in self.B])
        self.parspace_Atau = np.asarray(
          [(_A,_tau) for _A in self.A for _tau in self.tau])
        self.parspace_AB = np.asarray(
          [(_A,_B) for _A in self.A for _B in self.B])

        #Initialize quantities to be compute.
        self.signal = None

        #Used to compute the default tomography analysis.
        if case == 'default':
            self.subdir = 'test1/'    
            self.time_step = 2.4 #in units of second.
            #self.region_to_fit = 'all'
            self.region_to_fit = 'norest'
            self.show_fig = False
            self.save_fig = True

        if case == 'patient':
            self.subdir = 'patient/'    
            self.time_step = 2.4 #in units of second.
            #self.region_to_fit = 'all'
            self.region_to_fit = 'norest'
            self.show_fig = False
            self.save_fig = True
                
        create_run_dir('./../OUTPUT_FILES/RUNS/' + self.subdir)
        print '    N parspace cells = ', int(self.N_cells)
        print '    Done.\n'
            
