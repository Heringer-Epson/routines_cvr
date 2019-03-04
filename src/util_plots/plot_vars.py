#!/usr/bin/env python

import sys, os, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
from astropy import units as u

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'  
c = ['#fdae61', '#3288bd']
fs = 24.   

class Plot_Vars(object):
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
        self.fig, self.ax = plt.subplots(3,3, figsize=(16,16))
      
        self.tau, self.A, self.B, self.C, self.chi2 = None, None, None, None, None
        self.tau_est, self.A_est, self.B_est = None, None, None
        self.make_plot()
                
    def set_fig_frame(self):

        plt.subplots_adjust(
          left=.1, right=.95, bottom=.1, top=0.95, wspace=.15, hspace=0.15)
      
        self.ax[0,0].set_ylabel(r'$A\ \rm{[BOLD\ \%\ /\ mmHg]}$', fontsize=fs)
        self.ax[0,0].set_xlim(-0.05,0.05)
        self.ax[0,0].set_ylim(self._run.A_lim[0],self._run.A_lim[1])
        self.ax[0,0].tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax[0,0].tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax[0,0].tick_params('both', length=12, width=2., which='major',
                             right=True, top=True)
        self.ax[0,0].tick_params('both', length=6, width=2., which='minor',
                             right=True, top=True) 
        self.ax[0,0].xaxis.set_minor_locator(MultipleLocator(.01))
        self.ax[0,0].xaxis.set_major_locator(MultipleLocator(.05))
        self.ax[0,0].yaxis.set_minor_locator(MultipleLocator(.2))
        self.ax[0,0].yaxis.set_major_locator(MultipleLocator(1.))  
        self.ax[0,0].tick_params(labelbottom=False)  

        self.ax[1,0].set_ylabel(r'$\rm{log}\ \tau\ \rm{[s]}$', fontsize=fs)
        self.ax[1,0].set_xlim(-0.05,0.05)
        self.ax[1,0].set_ylim(self._run.logtau_lim[0],self._run.logtau_lim[1])
        self.ax[1,0].tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax[1,0].tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax[1,0].tick_params('both', length=12, width=2., which='major',
                             right=True, top=True)
        self.ax[1,0].tick_params('both', length=6, width=2., which='minor',
                             right=True, top=True) 
        self.ax[1,0].xaxis.set_minor_locator(MultipleLocator(.01))
        self.ax[1,0].xaxis.set_major_locator(MultipleLocator(.05))
        self.ax[1,0].yaxis.set_minor_locator(MultipleLocator(.2))
        self.ax[1,0].yaxis.set_major_locator(MultipleLocator(1.))  
        self.ax[1,0].tick_params(labelbottom=False)  

        self.ax[2,0].set_xlabel(r'$C\ \rm{[BOLD\ \%\ /\ s]}$', fontsize=fs)        
        self.ax[2,0].set_ylabel(r'$B\ \rm{[BOLD\ \%]}$', fontsize=fs)
        self.ax[2,0].set_xlim(-0.05,0.05)
        self.ax[2,0].set_ylim(self._run.B_lim[0],self._run.B_lim[1])
        self.ax[2,0].tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax[2,0].tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax[2,0].tick_params('both', length=12, width=2., which='major',
                             right=True, top=True)
        self.ax[2,0].tick_params('both', length=6, width=2., which='minor',
                             right=True, top=True) 
        self.ax[2,0].xaxis.set_minor_locator(MultipleLocator(.01))
        self.ax[2,0].xaxis.set_major_locator(MultipleLocator(.05))
        self.ax[2,0].yaxis.set_minor_locator(MultipleLocator(10.))
        self.ax[2,0].yaxis.set_major_locator(MultipleLocator(50.)) 

        self.ax[1,1].set_xlim(self._run.A_lim[0],self._run.A_lim[1])
        self.ax[1,1].set_ylim(self._run.logtau_lim[0],self._run.logtau_lim[1])
        self.ax[1,1].tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax[1,1].tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax[1,1].tick_params('both', length=12, width=2., which='major',
                             right=True, top=True)
        self.ax[1,1].tick_params('both', length=6, width=2., which='minor',
                             right=True, top=True) 
        self.ax[1,1].xaxis.set_minor_locator(MultipleLocator(.2))
        self.ax[1,1].xaxis.set_major_locator(MultipleLocator(1.))
        self.ax[1,1].yaxis.set_minor_locator(MultipleLocator(.2))
        self.ax[1,1].yaxis.set_major_locator(MultipleLocator(1.))  
        self.ax[1,1].tick_params(labelbottom=False, labelleft=False)  

        self.ax[2,1].set_xlabel(r'$A\ \rm{[BOLD\ \%\ /\ mmHg]}$', fontsize=fs)
        self.ax[2,1].set_xlim(self._run.A_lim[0],self._run.A_lim[1])
        self.ax[2,1].set_ylim(self._run.B_lim[0],self._run.B_lim[1])
        self.ax[2,1].tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax[2,1].tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax[2,1].tick_params('both', length=12, width=2., which='major',
                             right=True, top=True)
        self.ax[2,1].tick_params('both', length=6, width=2., which='minor',
                             right=True, top=True) 
        self.ax[2,1].xaxis.set_minor_locator(MultipleLocator(.2))
        self.ax[2,1].xaxis.set_major_locator(MultipleLocator(1.))
        self.ax[2,1].yaxis.set_minor_locator(MultipleLocator(10.))
        self.ax[2,1].yaxis.set_major_locator(MultipleLocator(50.))  
        self.ax[2,1].tick_params(labelleft=False)  

        self.ax[2,2].set_xlabel(r'$\rm{log}\ \tau\ \rm{[s]}$', fontsize=fs)
        self.ax[2,2].set_xlim(self._run.logtau_lim[0],self._run.logtau_lim[1])
        self.ax[2,2].set_ylim(self._run.B_lim[0],self._run.B_lim[1])
        self.ax[2,2].tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax[2,2].tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax[2,2].tick_params('both', length=12, width=2., which='major',
                             right=True, top=True)
        self.ax[2,2].tick_params('both', length=6, width=2., which='minor',
                             right=True, top=True) 
        self.ax[2,2].xaxis.set_minor_locator(MultipleLocator(.2))
        self.ax[2,2].xaxis.set_major_locator(MultipleLocator(1.))
        self.ax[2,2].yaxis.set_minor_locator(MultipleLocator(10.))
        self.ax[2,2].yaxis.set_major_locator(MultipleLocator(50.))  
        self.ax[2,2].tick_params(labelleft=False)  

        self.fig.delaxes(self.ax[0,1])
        self.fig.delaxes(self.ax[0,2])
        self.fig.delaxes(self.ax[1,2])
        
    def retrieve_data(self):
        fpath = ('./../OUTPUT_FILES/RUNS/' + self._run.subdir
                 + 'most_likely_pars.csv')
        self.M = pd.read_csv(fpath, header=0, low_memory=False)

        fpath = ('./../OUTPUT_FILES/RUNS/' + self._run.subdir
                 + 'estimated_A_tau_B.csv')
        self.E = pd.read_csv(fpath, header=0, low_memory=False)

    def plot_models(self):
        self.ax[0,0].errorbar(
          self.M['C'].values, self.M['A'].values,
          ls='None', marker='o', markersize=.2)

        self.ax[1,0].errorbar(
          self.M['C'].values, np.log10(self.M['tau'].values), 
          ls='None', marker='o', markersize=.2)

        self.ax[2,0].errorbar(
          self.M['C'].values, self.M['B'].values, 
          ls='None', marker='o', markersize=.2)

        self.ax[1,1].errorbar(
          self.M['A'].values, np.log10(self.M['tau'].values),
          ls='None', marker='o', markersize=.2)
          
        self.ax[2,1].errorbar(
          self.M['A'].values, self.M['B'].values,
          ls='None', marker='o', markersize=.2)

        self.ax[2,2].errorbar(
          np.log10(self.M['tau'].values), self.M['B'].values,
          ls='None', marker='o', markersize=.2)

    def manage_output(self):
        if self._run.save_fig:
            fpath = os.path.join(
              './../OUTPUT_FILES/RUNS/' + self._run.subdir, 'FIGURES/Fig_vars.pdf')
            plt.savefig(fpath, format='pdf')
        if self._run.show_fig:
            plt.show()
        plt.close()
            
    def make_plot(self):
        self.set_fig_frame()
        self.retrieve_data()
        self.plot_models()
        self.manage_output()

