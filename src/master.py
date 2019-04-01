#!/usr/bin/env python

from input_params import Input_Parameters
from define_domain import Define_Domain
from collect_data import Collect_Data
from norm_and_smooth import Norm_Smooth
from make_tau_interpolator import Tau_Interpolator
from estimate_bestpars import Estimate_Bestpars
from run_MCMC import Compute_Likelihoods
from store_data import Store_Data
from main_plotter import Main_Plotter

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
    
    def __init__(
      self, case, dirpath, rest, norest, step, ramp, do_inspect_domain, 
      do_collect_data, do_smooth, do_run_models, do_compute_likelihood,
      do_make_plots):
        
        self.do_inspect_domain = do_inspect_domain
        self.do_collect_data = do_collect_data
        self.do_smooth = do_smooth
        self.do_run_models = do_run_models
        self.do_compute_likelihood = do_compute_likelihood
        self.do_make_plots = do_make_plots
        
        self.run = Input_Parameters(case,dirpath,rest,norest,step,ramp)
        self.run_master()

    def run_master(self):       

        if self.do_collect_data:
            Collect_Data(self.run)
            if self.do_inspect_domain:
                Define_Domain(self.run)
        if self.do_smooth:
            Norm_Smooth(self.run)
        if self.do_run_models:
            Tau_Interpolator(self.run)
        if self.do_compute_likelihood:
            #Estimate_Bestpars(self.run)
            #Compute_Likelihoods(self.run)
            Store_Data(self.run)
        if self.do_make_plots:
            Main_Plotter(self.run)            
                    
if __name__ == '__main__':  
    Master(
      case='c02', dirpath='./../data_test/c02/', rest=(0., 110.),
      norest=(46.,390.), step=(46., 130.), ramp=(220.,390.), 
      do_inspect_domain=False, do_collect_data=False, do_smooth=False,
      do_run_models=False, do_compute_likelihood=False, do_make_plots=True)

    #Master(
    #  case='default', dirpath='./../data_test/normal1/', rest=(0., 110.),
    #  norest=(70.,730.), step=(150., 240.), ramp=(460.,670.), 
    #  do_inspect_domain=False, do_collect_data=True, do_smooth=False,
    #  do_run_models=False, do_compute_likelihood=False, do_make_plots=False)

    #Master(
    #  case='patient', dirpath='./../data_test/patient1/', rest=(0., 110.),
    #  norest=(70.,730.), step=(150., 240.), ramp=(460.,670.), 
    #  do_inspect_domain=False, do_collect_data=False, do_smooth=False,
    #  do_run_models=True, do_compute_likelihood=True, do_make_plots=False)

    #Master(
    #  case='tau', dirpath='./../data_test/tau-patient/', rest=(0., 120.),
    #  norest=(70.,730.), step=(150., 240.), ramp=(460.,670.), 
    #  do_inspect_domain=False, do_collect_data=False, do_smooth=False,
    #  do_run_models=False, do_compute_likelihood=False, do_make_plots=True)
