#!/usr/bin/env python
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import cPickle
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
import data_handling 


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
        print 'Plotting signal and best model...'
                
        self._run = _run
        self._idx_space = _idx_space
        self.S, self.I = None, None
        self.A_best, self.tau_best, self.B_best = None, None, None 
        self.A_est, self.tau_est, self.B_est = None, None, None 
        
        self.fig = plt.figure(figsize=(10,10))
        self.ax = self.fig.add_subplot(111)
                
        self.make_plot()
        
    def set_fig_frame(self):
        x_label = r'$t\ (\rm{s})$'
        y_label = r'BOLD$\ (\%)$'
        self.ax.set_xlabel(x_label, fontsize=20.)
        self.ax.set_ylabel(y_label, fontsize=20.)
        self.ax.set_xlim(0., 850.)
        self.ax.set_ylim(80., 130.)
        self.ax.tick_params(axis='y', which='major', labelsize=20., pad=8)      
        self.ax.tick_params(axis='x', which='major', labelsize=20., pad=8)
        self.ax.minorticks_off()
        self.ax.tick_params(
          'both', length=8, width=1., which='major', direction='in')
        self.ax.tick_params(
          'both', length=4, width=1., which='minor', direction='in')    
        self.ax.xaxis.set_minor_locator(MultipleLocator(50.))
        self.ax.xaxis.set_major_locator(MultipleLocator(100.))
        self.ax.yaxis.set_minor_locator(MultipleLocator(1.))
        self.ax.yaxis.set_major_locator(MultipleLocator(5.))  

    def retrieve_data(self):

        #Load observational data.
        fpath_S = './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'PICKLES/smooth.pkl'
        fS = open(fpath_S, 'r')
        self.S = cPickle.load(fS)
        fS.close()        

        #Load interpolator for computing pCO2 convolutions.
        fpath_I = './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'PICKLES/tau_interp.pkl'
        fI = open(fpath_I, 'r')
        self.I = cPickle.load(fI)
        fI.close()
        
        #Load best pars.
        fpath = './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'most_likely_pars.csv'
        tau, A, B = np.loadtxt(
          fpath, skiprows=1, delimiter=',', usecols=(1,4,7), unpack=True)
        self.tau_best, self.A_best, self.B_best =\
          tau[self._idx_space], A[self._idx_space], B[self._idx_space]
          #tau[self._idx_space], A[self._idx_space], B[self._idx_space]

        fpath_est = './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'estimated_A_tau_B.csv'
        tau, A, B = np.loadtxt(
          fpath_est, skiprows=1, delimiter=',', usecols=(1,2,3), unpack=True)
        self.tau_est, self.A_est, self.B_est =\
          tau[self._idx_space], A[self._idx_space], B[self._idx_space]

    def plot_quantities(self):

        t, pCO2 = self.S['time'], self.S['pCO2']
        y = self.S['signal_ns'][self._idx_space,:]
        yerr = self.S['signal_noise'][self._idx_space,:]
        pCO2 = self.S['pCO2']
        
        pCO2_conv_est = self.I(t, self.tau_est)
        y_model_est = self.A_est * pCO2_conv_est + self.B_est

        pCO2_conv_best = self.I(t, self.tau_best)
        y_model_best = self.A_best * pCO2_conv_best + self.B_best

        self.ax.errorbar(
          t, y, yerr=yerr, ls='None', marker='o', color='k',
          alpha=0.6, markersize=5., capsize=0., elinewidth=0.4, label=r'Signal')          
        self.ax.plot(t, y_model_est, ls='--', marker='None', color='darkgreen',
                     lw=2., label=r'Est fit')   
        self.ax.plot(t, y_model_best, ls='--', marker='None', color='firebrick',
                     lw=2., label=r'Best fit')   

        #self.ax.plot(t, pCO2 + self.B_est, ls='-', marker='None',  lw=2., color='m',
        #             label=r'pCO2$\ (\rm{mmHg})$')
               
        self.ax.legend(frameon=False, fontsize=20., numpoints=1, loc='best')            

    def manage_output(self):
        if self._run.save_fig:
            fpath = os.path.join('./../OUTPUT_FILES/RUNS/' + self._run.subdir,
                                 'FIGURES/Fig_most_likely_A_tau.pdf')
            plt.savefig(fpath, format='pdf')
        if self._run.show_fig:
            plt.show()
        plt.close()

    def make_plot(self):
        self.set_fig_frame()
        self.retrieve_data()
        self.plot_quantities()
        self.manage_output()
        print '    Done.\n'
