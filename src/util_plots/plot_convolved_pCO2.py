#!/usr/bin/env python

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import cPickle
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
import stats
import data_handling 

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
fs = 24

class Plot_Convolved(object):
    """
    Description:
    ------------
    Plot a series of convolved pCO2 (with an exponential) for a series of taus.
    This is a visual test of whether the interpolator code
    (./../make_tau_interpolator.py) is properly working or not.
    
    Parameters:
    -----------
    _run : object
      instance of input_params containing the data in each case of study.
    
    Outputs:
    --------
    ('./../OUTPUT_FILES/RUNS/' + self._run.subdir + 'Figures/'
      + 'Fig_convolved_pCO2.pdf')
    """        
    def __init__(self, _run):
        print 'Plotting confidence contour...'
        self._run = _run

        fpath_D = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'PICKLES/data.pkl'
        fpath_I = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'PICKLES/tau_interp.pkl'
        fpath_S = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'PICKLES/smooth.pkl'
        fD, fS, fI = open(fpath_D , 'r'), open(fpath_S , 'r'), open(fpath_I , 'r')
        self.D, self.S, self.I = cPickle.load(fD), cPickle.load(fS), cPickle.load(fI) 
        fD.close(), fS.close(), fI.close()
        
        self.fig = plt.figure(figsize=(10,10))
        self.ax = self.fig.add_subplot(111)
        
        self.make_plot()
        
    def set_fig_frame(self):
        x_label = r'$t\, [\rm{s}]$'
        y_label = r'$[\rm{mmHg}]$'
        self.ax.set_xlabel(x_label, fontsize=fs)
        self.ax.set_ylabel(y_label, fontsize=fs)
        self.ax.set_xlim(-550., 900.)
        self.ax.set_ylim(30., 55.)
        self.ax.tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax.tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax.minorticks_off()
        self.ax.tick_params(
          'both', length=8, width=1., which='major', direction='in')
        self.ax.tick_params(
          'both', length=4, width=1., which='minor', direction='in')    
        self.ax.xaxis.set_minor_locator(MultipleLocator(50))
        self.ax.xaxis.set_major_locator(MultipleLocator(250))
        self.ax.yaxis.set_minor_locator(MultipleLocator(1.))
        self.ax.yaxis.set_major_locator(MultipleLocator(10.))
    
    def plot_convolved_pCO2(self):
        t, pCO2 = self.D['time'], self.D['pCO2']
        t_aux = np.arange(-500., 800, 2.4)

        self.ax.errorbar(
          self.S['time'], self.S['pCO2'], yerr=self.S['pCO2_noise'], ls='None',
          marker='o', color='m', alpha=0.6, markersize=2., capsize=0.,
          elinewidth=0.4, zorder=3) 

        for tau in self._run.tau:
            self.ax.plot(
              t_aux,self.I(t_aux, tau),marker='None',ls='-',color='r',
              lw=2., alpha=0.7, zorder=2)  

        for tau in np.logspace(0., 2., 100):
            self.ax.plot(t_aux,self.I(t_aux, tau),marker='None',ls='-',color='grey',
                        lw=1., alpha=0.2, zorder=1) 

    def make_legend(self):
        self.ax.plot(
          [np.nan], [np.nan], marker='None',ls='-',color='grey', lw=1.,
          alpha=1., label=r'Interpolated Models')
        self.ax.plot(
          [np.nan], [np.nan], marker='None',ls='-',color='r', lw=2.,
          label=r'Convolved $\rm{pCO_2}$')
        self.ax.errorbar(
          [np.nan], [np.nan], yerr=[np.nan], ls='None',
          marker='o', color='m', alpha=1., markersize=5., capsize=0.,
          elinewidth=1., label=r'$\rm{pCO_2}$ Reading') 

        self.ax.legend(frameon=False, fontsize=fs, numpoints=1, loc=2)            
        

    def manage_output(self):
        if self._run.save_fig:
            fpath = os.path.join('./../OUTPUT_FILES/RUNS/' + self._run.subdir,
                                 'FIGURES/Fig_convolved_pCO2.pdf')
            plt.savefig(fpath, format='pdf')
        if self._run.show_fig:
            plt.show() 
        plt.close()
        
    def make_plot(self):
        self.set_fig_frame()
        self.plot_convolved_pCO2()
        self.make_legend()
        self.manage_output()
