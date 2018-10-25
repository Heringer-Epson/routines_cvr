#!/usr/bin/env python

from util_plots.plot_models import Plot_Models


class Main_Plotter(object):
    """
    Description:
    ------------
    This piece of code works as an interface between master.py and modules
    that produce relevant plots for each each individual run.

    Parameters:
    -----------
    _inputs : ~instance
        Instance of the Input_Parameters class defined in input_params.py.
    """    
    def __init__(self, _run):
        pass
        #Plot_Models(_run)


