#!/usr/bin/env python

import os, sys
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'

class Plot_Signal(object):
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
    def __init__(self, _time, _signal, _model, _pCO2):
                
        self._time = _time
        self._signal = _signal
        self._model = _model
        self._pCO2 = _pCO2

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
        self.ax.plot(self._time, self._signal, ls='-', marker='None', color='k', 
                     label=r'Signal')
        self.ax.plot(self._time, self._pCO2, ls='-', marker='None', color='m', 
                     label=r'pCO2$\ (\rm{mmHg})$')
        self.ax.plot(self._time, self._model, ls='--', marker='None', color='firebrick', 
                     label=r'Best fit')      
        self.ax.legend(frameon=True, fontsize=20., numpoints=1, loc='best')            

    def make_plot(self):
        self.set_fig_frame()
        self.plot_quantities()
        plt.show()  
