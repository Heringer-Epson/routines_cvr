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

        self.A_lim = (-1.,1.)
        self.B_lim = (20.,180.)
        self.C_lim = (-100.,100.)
        self.tau_lim = (1.,500.)
        
        #self.A_lim = (-2.5,2.5)
        #self.B_lim = (-50.,200.)
        #self.C_lim = (-1.,1.)
        #self.tau_lim = (0.1,1000.)
        self.logtau_lim = (np.log10(self.tau_lim[0]),np.log10(self.tau_lim[1]))
        self.tau = np.logspace(self.logtau_lim[0], self.logtau_lim[1], 100)
        self.t_pivot = 400.

        #Initialize quantities to be compute.
        self.signal = None

        #Used to compute the default tomography analysis.
        if case == 'default':
            self.subdir = 'healthy/'    
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

        if case == 'tau':
            self.subdir = 'tau-patient/'    
            self.time_step = 2.4 #in units of second.
            #self.region_to_fit = 'all'
            self.region_to_fit = 'norest'
            self.show_fig = False
            self.save_fig = True

        if case == 'c02':
            self.subdir = 'c02/'    
            self.time_step = 2.4 #in units of second.
            #self.region_to_fit = 'all'
            self.region_to_fit = 'norest'
            self.show_fig = False
            self.save_fig = True
                
        create_run_dir('./../OUTPUT_FILES/RUNS/' + self.subdir)
        print '    Done.\n'
            
