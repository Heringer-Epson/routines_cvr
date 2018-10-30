#!/usr/bin/env python

import os, sys
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
from lib import data_handling 


mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

class Plot_Signal(object):
    """
    Description:
    ------------
    Makes a figure showing, for a given voxel, the input pCO2 stimulus, the
    output signal and the best model.
    
    Parameters
    ----------
    _run : object
      instance of input_params containing the data in each case of study.
    _idx_space : ~int
      Spatial index of the voxel to be plotted.
    
    Outputs:
    --------
    ('./../OUTPUT_FILES/RUNS/' + self._run.subdir + 'Figures/'
      + 'Fig_most_likely_A_tau_B.csv')
    """        
    def __init__(self, _run, _idx_space):
                
        self._run = _run
        self._idx_space = _idx_space

        self.fig = plt.figure(figsize=(10,10))
        self.ax = self.fig.add_subplot(111)
        
        self._model = None
        
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

        set_models = data_handling.get_convolved_models(self._run)
        
        fpath = os.path.join(
          './../OUTPUT_FILES/RUNS/' + self._run.subdir, 'most_likely_A_tau_B.csv')

        df = pd.read_csv(fpath, header=0, low_memory=False, dtype='str')            
        cond = (df['voxel'] == str(self._idx_space)).values        
        
        A = float(df['A'].values[cond][0])
        tau = float(df['tau'].values[cond][0])
        B = float(df['B'].values[cond][0])
        self.model = set_models[data_handling.pars2label(A,tau,B)]

    def plot_quantities(self):
        self.ax.plot(self._run.time[1:], self._run.signal_ns[3523,:][1:], ls='-',
                     marker='None', color='k', label=r'Signal')
        self.ax.plot(self._run.time[1:], self._run.pCO2[1:], ls='-',
                     marker='None', color='m', label=r'pCO2$\ (\rm{mmHg})$')
        self.ax.plot(self._run.time[1:], self.model[1:], ls='--', marker='None',
                     color='firebrick', label=r'Best fit')      
        self.ax.legend(frameon=True, fontsize=20., numpoints=1, loc='best')            

    def manage_output(self):
        if self._run.save_fig:
            fpath = os.path.join('./../OUTPUT_FILES/RUNS/' + self._run.subdir,
                                 'Figures/Fig_most_likely_A_tau_B.csv')
            plt.savefig(fpath, format='pdf')
        if self._run.show_fig:
            plt.show() 

    def make_plot(self):
        self.set_fig_frame()
        self.get_most_likely_model()
        self.plot_quantities()
        self.manage_output()
