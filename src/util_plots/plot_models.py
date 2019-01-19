#!/usr/bin/env python
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
import data_handling 

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

class Plot_Models(object):
    """
    Description:
    ------------
    TBW
    Parameters:
    -----------
    TBW
    
    Outputs:
    --------
    None
    """        
    def __init__(self, _run, _idx_space):
                
        self._run = _run
        self._idx_space = _idx_space
        self.M = data_handling.get_master_pickle(self._run)        

        self.fig = plt.figure(figsize=(10,10))
        self.ax = self.fig.add_subplot(111)
        
        self.most_likely = None
        
        self.make_plot()
        
    def set_fig_frame(self):
        x_label = r'$t\ (\rm{s})$'
        y_label = r'BOLD$\ (\%)$'
        self.ax.set_xlabel(x_label, fontsize=20.)
        self.ax.set_ylabel(y_label, fontsize=20.)
        #self.ax.set_xlim(14., 18.)
        #self.ax.set_ylim(0., 1.6)
        self.ax.tick_params(axis='y', which='major', labelsize=20., pad=8)      
        self.ax.tick_params(axis='x', which='major', labelsize=20., pad=8)
        self.ax.minorticks_off()
        self.ax.tick_params(
          'both', length=8, width=1., which='major', direction='in')
        self.ax.tick_params(
          'both', length=4, width=1., which='minor', direction='in')    
        #self.ax.xaxis.set_minor_locator(MultipleLocator(.5))
        #self.ax.xaxis.set_major_locator(MultipleLocator(1.))
        #self.ax.yaxis.set_minor_locator(MultipleLocator(.05))
        #self.ax.yaxis.set_major_locator(MultipleLocator(.2))  

    def get_most_likely_model(self):
        
        fpath = os.path.join(
          './../OUTPUT_FILES/RUNS/' + self._run.subdir, 'most_likely_A_tau.csv')

        df = pd.read_csv(fpath, header=0, low_memory=False, dtype='str')            
        cond = (df['voxel'] == str(self._idx_space)).values        
                
        A = float(df['A'].values[cond][0])
        tau = float(df['tau'].values[cond][0])
        self.most_likely = self.M[data_handling.pars2label(A,tau)]
    
    def plot_quantities(self):

        self.ax.errorbar(
          self.M['time'][1:], self.M['signal_der'][3523,:],
          yerr=self.M['signal_noise'][3523,:], ls='None', marker='o', color='k',
          alpha=0.6, label=r'Signal')        
        
        self.ax.plot(self.M['time'][1:], self.M['signal_der_s'][3523,:], ls='-',
                     marker='None', color='k', label=r'Signal')
        
        #self.ax.plot(self.M['time'][1:], self.most_likely, ls='--', marker='None',
        #             color='firebrick', label=r'Best fit')   
        
                     
        #self.ax.plot(
        #  self._run.time, self._run.pCO2, ls='-', marker='None', color='m', 
        #  alpha=0.5, label=r'pCO2$\ (\rm{mmHg})$')
        #self.ax.errorbar(
        #  self._run.time, self._run.pCO2, yerr=self._run.pCO2_noise, ls='None',
        #  marker='o', markersize=2., color='m', elinewidth=0.5, capsize=0.,
        #  label=r'Signal (smoothed)')

    def add_models(self):
        models = data_handling.get_master_pickle(self._run)

        for (A,tau) in self.M['parspace']:
            model = self.M[data_handling.pars2label(A,tau)]
            self.ax.plot(self.M['time'][1:], model, ls='--', marker='None',
            color='gray', alpha=0.3)      
       
        #self.ax.legend(frameon=True, fontsize=20., numpoints=1, loc='best')            

    def manage_output(self):
        if self._run.save_fig:
            fpath = os.path.join('./../OUTPUT_FILES/RUNS/' + self._run.subdir,
                                 'FIGURES/Fig_models.pdf')
            plt.savefig(fpath, format='pdf')
        if self._run.show_fig:
            plt.show()
        plt.close()

    def make_plot(self):
        self.set_fig_frame()
        self.get_most_likely_model()
        self.plot_quantities()
        self.add_models()
        self.manage_output()
