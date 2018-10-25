#!/usr/bin/env python

import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'src'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
from lib import stats

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

class Plot_Contour(object):
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
    def __init__(self, _parspace, _ln_L):
                
        self._parspace = _parspace
        self._ln_L = _ln_L

        self.fig = plt.figure(figsize=(10,10))
        self.ax = self.fig.add_subplot(111)
        
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
    
    def plot_quantities(self):
        columns = zip(*self._parspace)
        x, y = np.array(columns[0]), np.array(columns[1])
        #print len(x), len(y)
        stats.plot_contour(self.ax, x, y, self._ln_L, 'b', add_max=True)    
        
        #self.ax.legend(frameon=True, fontsize=20., numpoints=1, loc='best')            

    def make_plot(self):
        #self.set_fig_frame()
        self.plot_quantities()
        plt.show()  
