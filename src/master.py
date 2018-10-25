#!/usr/bin/env python

import os
from input_params import Input_Parameters
from get_data import Read_Data
from norm_and_smooth import Norm_Smooth
from define_domain import Define_Domain
from determine_baseline import Fit_Baseline
from make_A_tau_model import Save_Models
from compute_likelihood import Compute_Likelihoods

class Master(object):
    """
    Code Description
    ----------    
    TBW.

    Parameters
    ----------
    TBW.
    """
    
    def __init__(self, case, dirpath, rest, step, ramp, inspect_domain, 
                 run_models, show_fig, save_fig):
        
        #Flags.
        self.inspect_domain = inspect_domain
        self.run_models = run_models
        
        #Initialize run (model object which stores relevant quantities).
        self.run = Input_Parameters(case,dirpath,rest,step,ramp,show_fig,save_fig)
            
        self.run_master()

    def run_master(self):       
        
        self.run.signal, self.run.pCO2, self.run.time =\
          Read_Data(self.run).run_task()
        self.run.signal_n, self.run.signal_ns, self.run.signal_sn_unc,\
          self.run.pCO2_s, self.run.pCO2_noise = Norm_Smooth(self.run).run_task()

        if self.inspect_domain:
            Define_Domain(self.run)
        else:
            if self.run_models:
                Save_Models(self.run)                                          
            #self.run.B, self.run.B_unc = Fit_Baseline(self.run).run_task()
            Compute_Likelihoods(self.run, 'step').run_task()
            
                    
if __name__ == '__main__':  
    Master('default', './../data_test/normal1/', (0., 110.), (150., 240.), 
           (460.,670.), False, False, True, False)


