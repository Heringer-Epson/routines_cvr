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

class Plot_Contour(object):
    """
    Description:
    ------------
    Plot the likelihood contour of a given voxel.
    
    Parameters:
    -----------
    _run : object
      instance of input_params containing the data in each case of study.
    _idx_space : ~int
      Spatial index of the voxel to be plotted.
    
    Outputs:
    --------
    ('./../OUTPUT_FILES/RUNS/' + self._run.subdir + 'Figures/'
      + 'Fig_A_tau_contour.csv')
    """        
    def __init__(self, _run, _idx_space):
        print 'Plotting confidence contour...'
        self._run = _run
        self._idx_space = _idx_space

        fpath_L = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'PICKLES/likelihoods.pkl'
        fpath_S = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'PICKLES/smooth.pkl'
        fL, fS = open(fpath_L , 'r'), open(fpath_S , 'r')
        self.L, self.S = cPickle.load(fL), cPickle.load(fS) 
        fL.close(), fS.close()
        
        self.fig = plt.figure(figsize=(10,10))
        self.ax = self.fig.add_subplot(111)
        
        self.make_plot()
        
    def set_fig_frame(self):
        x_label = r'$ssCVR$'
        y_label = r'$\tau\, [\rm{s}]$'
        self.ax.set_yscale('log')
        self.ax.set_xlabel(x_label, fontsize=20.)
        self.ax.set_ylabel(y_label, fontsize=20.)
        self.ax.set_xlim(-1., 1.)
        self.ax.set_ylim(0.01, 100.)
        self.ax.tick_params(axis='y', which='major', labelsize=20., pad=8)      
        self.ax.tick_params(axis='x', which='major', labelsize=20., pad=8)
        self.ax.minorticks_off()
        self.ax.tick_params(
          'both', length=8, width=1., which='major', direction='in')
        self.ax.tick_params(
          'both', length=4, width=1., which='minor', direction='in')    
        self.ax.xaxis.set_minor_locator(MultipleLocator(.1))
        self.ax.xaxis.set_major_locator(MultipleLocator(.5))
        locmaj = mpl.ticker.LogLocator(base=10,numticks=12) 
        locmin = mpl.ticker.LogLocator(
          base=10.0,subs=np.arange(0.1,0.91,0.1),numticks=12)
        self.ax.yaxis.set_major_locator(locmaj)
        self.ax.yaxis.set_minor_locator(locmin)
        self.ax.yaxis.set_minor_formatter(mpl.ticker.NullFormatter())
    
    def plot_quantities(self):
        ln_L = self.L['likelihood_list_' + str(self._idx_space)]
        columns = zip(*self._run.parspace_Atau)
        x, y = np.array(columns[0]), np.array(columns[1])
        stats.plot_contour(self.ax, x, y, ln_L, 'b', add_max=True)    

    def manage_output(self):
        if self._run.save_fig:
            fpath = os.path.join('./../OUTPUT_FILES/RUNS/' + self._run.subdir,
                                 'FIGURES/Fig_A_tau_contour.csv')
            plt.savefig(fpath, format='pdf')
        if self._run.show_fig:
            plt.show() 
        plt.close()
        
    def make_plot(self):
        self.set_fig_frame()
        self.plot_quantities()
        self.manage_output()
