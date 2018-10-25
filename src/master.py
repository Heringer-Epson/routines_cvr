#!/usr/bin/env python

from input_params import Input_Parameters
from norm_and_smooth import Norm_Smooth
from define_domain import Define_Domain
from make_A_tau_model import Save_Models
from compute_likelihood import Compute_Likelihoods
from main_plotter import Main_Plotter

from lib.data_handling import Read_Data

class Master(object):
    """
    Code Description
    ----------    
    Main code. It call all major routines to read, the input data, create A-tau
    models, find the most likely parameters and make relevant figure. The
    information produced is stored as a RUN.

    Parameters
    ----------
    case : ~str
        A string keyword to be used to retrieve a set of input parameters from
        Input_Parameters. E.g. "case='default'".
    dirpath : ~str
        Relative path to the directory containing the patient data to be
        analysed. E.g. "dirpath='./../data_test/normal1/'".
    rest : ~tuple
        Tuple with two floats, containing the starting and ending time of the
        "rest" phase of the pCO2 stimulus. This can be determined in a pre-run,
        where the flag "inspect_domain" (see below) is set equal True (bool).
    step : ~tuple
        Same as "rest", but for the step phase.
    ramp : ~tuple
        Same as "rest", but for the ramp phase.    
    inspect_domain : ~bool
        If True, this code will show a figure of the pCO2 stimulus, so that
        the rest, step and ramp phases can be determined by eye. This flag is 
        recommended to be set True every time a new patient is studied. After
        these phases are fixed once, one may set the flag False henceforth.
    run_models : ~bool
        If True, a routine that creates A-tau models will be invoked and a set
        of models (according to parameters in input_params) will be created
        and stored in a .pkl file.
    plots_flag : ~bool
        If True, the main_plotter routine will be called and a set of figures
        will be created.

    Outputs:
    --------
    './../OUTPUT_FILES/RUNS/' + subdir (where subdir is set in input_params for
    each "case").

    """
    
    def __init__(self, case, dirpath, rest, step, ramp, inspect_domain, 
                 run_models, plots_flag):
        
        self.inspect_domain = inspect_domain
        self.run_models = run_models
        self.plots_flag = plots_flag
        
        self.run = Input_Parameters(case,dirpath,rest,step,ramp)
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
            Compute_Likelihoods(self.run, 'step').run_task()

        if self.plots_flag:
            Main_Plotter(self.run)            
                    
if __name__ == '__main__':  
    Master(case='default', dirpath='./../data_test/normal1/', rest=(0., 110.),
           step=(150., 240.), ramp=(460.,670.), inspect_domain=False,
           run_models=False, plots_flag=False)
