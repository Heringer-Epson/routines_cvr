#!/usr/bin/env python

import sys, os, time
import numpy as np
import pandas as pd
import corner
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
from astropy import units as u

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'  
c = ['#fdae61', '#3288bd']
fs = 14.  

class Plot_Corner(object):
    """
    Description:
    ------------
    TBW.

    Parameters:
    -----------
    show_fig : ~bool
        True of False. Whether to show or not the produced figure.
    save_fig : ~bool
        True of False. Whether to save or not the produced figure.

    Notes:
    ------
    The normalization used here is different than in Heringer+ 2017. In that
    paper the DTD is normalized at 0.5 Gyr, whereas here an arbitraty
    constant (A=10**-12) is given to the DTD of the form SN rate = A*t**s1
    where s1 is the slope prior to 1 Gyr. 
             
    Outputs:
    --------
    ./../OUTPUT_FILES/ANALYSES_FIGURES/Fig_sSNRL.pdf
    
    References:
    -----------
    Heringer+ 2017: http://adsabs.harvard.edu/abs/2017ApJ...834...15H
    """         
    def __init__(self, _run):
        self._run = _run
        self.fig, self.ax = plt.subplots(2,2, figsize=(12,12))
      
        self.tau, self.A, self.B, self.C, self.chi2 = None, None, None, None, None
        self.tau_est, self.A_est, self.B_est = None, None, None
        self.make_plot()
            
        
    def retrieve_data(self):
        fpath = ('./../OUTPUT_FILES/RUNS/' + self._run.subdir
                 + 'most_likely_pars.csv')
        self.M = pd.read_csv(fpath, header=0, low_memory=False)

        fpath = ('./../OUTPUT_FILES/RUNS/' + self._run.subdir
                 + 'estimated_A_tau_B.csv')
        self.E = pd.read_csv(fpath, header=0, low_memory=False)

    def make_corner_plot(self):
        
        data = np.column_stack((
          self.M['C'].values, self.M['A'].values, 
          np.log10(self.M['tau'].values), self.M['B'].values))

        figure = corner.corner(
          data, labels=[r'$C\ \rm{[BOLD\ \%\ /\ s]}$',
          r'$A\ \rm{[BOLD\ \%\ /\ mmHg]}$', r'$\rm{log}\ \tau\ \rm{[s]}$',
          r'$B\ \rm{[BOLD\ \%]}$'], quantiles=[0.16, 0.5, 0.84],
          show_titles=True, title_kwargs={'fontsize': fs},
          label_kwargs={'fontsize': fs}) 
        
    def manage_output(self):
        if self._run.save_fig:
            fpath = os.path.join(
              './../OUTPUT_FILES/RUNS/' + self._run.subdir, 'FIGURES/Fig_corner.pdf')
            plt.savefig(fpath, format='pdf')
        if self._run.show_fig:
            plt.show()
        plt.close()
            
    def make_plot(self):
        self.retrieve_data()
        self.make_corner_plot()
        self.manage_output()

