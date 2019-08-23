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

pars2file = {'est': 'estimated_pars.csv', 'proc': 'most_likely_pars.csv',
             'procgmwm': 'most_likely_pars_gmwm.csv'}

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
    def __init__(self, _run, mode, gmwm=''):
        self._run = _run
        self.mode, self.gmwm = mode, gmwm

        self.fig, self.ax = plt.subplots(2,2, figsize=(12,12))
        self.make_plot()
        
    def retrieve_data(self):
        fname = pars2file[self.mode + self.gmwm]
        fpath = './../OUTPUT_FILES/RUNS/' + self._run.subdir + fname
        self.df = pd.read_csv(fpath, header=0, low_memory=False)

    def make_corner_plot(self):
                
        if self.gmwm == '':
            data = np.column_stack((
              self.df['C'].values, self.df['A'].values, 
              np.log10(self.df['tau'].values), self.df['B'].values))

            figure = corner.corner(
              data, labels=[r'$C\ \rm{[BOLD\ \%\ /\ s]}$',
              r'$A\ \rm{[BOLD\ \%\ /\ mmHg]}$', r'$\rm{log}\ \tau\ \rm{[s]}$',
              r'$B\ \rm{[BOLD\ \%]}$'], quantiles=[0.16, 0.5, 0.84],
              show_titles=True, title_kwargs={'fontsize': fs},
              label_kwargs={'fontsize': fs}) 

        if self.gmwm == 'gmwm':
            data = np.column_stack((
              self.df['C'].values, self.df['A_gm'].values, self.df['A_wm'].values, 
              np.log10(self.df['tau_gm'].values), np.log10(self.df['tau_wm'].values),
               self.df['B'].values))

            figure = corner.corner(
              data, labels=[r'$C\ \rm{[BOLD\ \%\ /\ s]}$',
              r'$A_{gm}\ \rm{[BOLD\ \%\ /\ mmHg]}$',  r'$A_{wm}\ \rm{[BOLD\ \%\ /\ mmHg]}$',
              r'$\rm{log}\ \tau_{gm}\ \rm{[s]}$', r'$\rm{log}\ \tau_{wm}\ \rm{[s]}$',
              r'$B\ \rm{[BOLD\ \%]}$'], quantiles=[0.16, 0.5, 0.84],
              show_titles=True, title_kwargs={'fontsize': fs},
              label_kwargs={'fontsize': fs}) 
        
    def manage_output(self):
        if self._run.save_fig:
            fpath = os.path.join(
              './../OUTPUT_FILES/RUNS/' + self._run.subdir,
              'FIGURES/Fig_corner' + self.gmwm + '.png')
            plt.savefig(fpath, format='png')
        if self._run.show_fig:
            plt.show()
        plt.close()
            
    def make_plot(self):
        self.retrieve_data()
        self.make_corner_plot()
        self.manage_output()

