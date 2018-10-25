import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
fs = 24.

class Define_Domain(object):
    """
    Description:
    ------------
    This routine should be executed prior to the full analysis in master.py. It
    is used to the visually define (in time space) the step and ramp for each
    case of study.  
    
    Parameters:
    -----------
    _run : object
      instance of input_params containing the data in each case of study.
    
    Outputs:
    --------
    None
    """   
        
    def __init__(self, _run):
        self._run = _run
        
        fig = plt.figure(figsize=(10,10))
        self.ax = fig.add_subplot(111)
        
        self.make_plot()

    def set_fig_frame(self):        
        x_label = r'$t\, \rm{[s]}$'
        y_label = r'$P(\rm{CO_2})\, \rm{[mmHg]}$'
        self.ax.set_xlabel(x_label, fontsize=fs)
        self.ax.set_ylabel(y_label, fontsize=fs)
        self.ax.set_xlabel(x_label, fontsize=fs)
        self.ax.set_ylabel(y_label, fontsize=fs)
        
        self.ax.tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax.tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax.minorticks_off()
        self.ax.tick_params(
          'both', length=8, width=1., which='major', direction='in')
        self.ax.tick_params(
          'both', length=4, width=1., which='minor', direction='in')
        self.ax.xaxis.set_ticks_position('both')
        self.ax.yaxis.set_ticks_position('both')             

        self.ax.xaxis.set_minor_locator(MultipleLocator(50))
        self.ax.xaxis.set_major_locator(MultipleLocator(100))
        self.ax.yaxis.set_minor_locator(MultipleLocator(1))
        self.ax.yaxis.set_major_locator(MultipleLocator(5))

    def plot_data(self):
        self.ax.plot(self._run.time, self._run.pCO2_s, ls='-', lw=1.5, color='r')

        self.ax.errorbar(
          self._run.time, self._run.pCO2, yerr=self._run.pCO2_noise, marker='o',
          ls='None', markersize=1.8, color='r', elinewidth=0.5, capsize=0.,
          alpha=0.4)

    def manage_output(self):
        if self._run.save_fig:
            fpath = os.path.join('./../OUTPUT_FILES/RUNS/' + self._run.subdir
                                 + 'FIGURES/', 'Fig_pCO2.pdf')
            plt.savefig(fpath, format='pdf')
        if self._run.show_fig:
            plt.show()      
        
    def make_plot(self):
        self.set_fig_frame()
        self.plot_data()
        self.manage_output()
