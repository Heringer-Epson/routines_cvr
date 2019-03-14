#!/usr/bin/env python

from util_plots.plot_convolved_pCO2 import Plot_Convolved
from util_plots.plot_signal import Plot_Signal
from util_plots.plot_models import Plot_Models
from util_plots.plot_likelihood import Plot_Contour
from util_plots.plot_brain import Plot_Brain
from util_plots.plot_histograms import Plot_Hist
from util_plots.plot_vars import Plot_Vars
from util_plots.plot_corner import Plot_Corner

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
        #Plot_Convolved(_run)
        #Plot_Signal(_run, 0)
        #Plot_Contour(_run, 3523)
        #Plot_Models(_run, 3523)
        
        #Plot_Hist(_run)
        #Plot_Vars(_run)
        #Plot_Corner(_run)
        #Plot_Brain(_run, 'A', 'est')
        Plot_Brain(_run, 'A', 'proc')
        #Plot_Brain(_run, 'tau', 'est')
        #Plot_Brain(_run, 'tau', 'proc')
        #Plot_Brain(_run, 'B', 'est')
        #Plot_Brain(_run, 'B', 'proc')
        #Plot_Brain(_run, 'A_unc', 'proc')
        pass

