#!/usr/bin/env python

import sys, os, time
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator
from astropy import units as u
from scipy.optimize import curve_fit

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'  
c = ['#e41a1c','#377eb8','#4daf4a']
fs = 24.   

def add_curves(ax,qtty,xlim,N_bins,c):
    step = (xlim[1] - xlim[0]) / N_bins
    bins = np.arange(xlim[0], xlim[1] + 1.e-5 ,step)
    ax.hist(qtty,bins,align='mid',normed=True,histtype='step',color=c,lw=2.)
    
class Plot_Hist(object):
    """
    Description:
    ------------
    Makes the Fig. 1 of the DTD paper, displaying SN rate (per unit of
    luminosity) as a function of Dcolor. 4 panels are included, showing the
    impact of choosing different parameters, such as the IMF, the metallicity,
    the SFH and the time of onset of SNe Ia. Similar replicates Fig. 3
    in Heringer+ 2017.

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
        self.fig, self.ax = plt.subplots(1,2, figsize=(14,8))
      
        self.df_est, self.df_mcmc = None, None

        self.make_plot()
                
    def set_fig_frame(self):

        plt.subplots_adjust(
          left=.07, right=.95, bottom=.15, top=0.9, wspace=.15, hspace=0.4)

        self.ax[0].text(
          -0.08, .5, r'Brain pdf', fontsize=fs, transform=self.ax[0].transAxes,
          rotation=90., horizontalalignment='center',
          verticalalignment='center')        
        self.ax[0].set_xlabel(r'$A\ \rm{[BOLD\ \%\ /\ mmHg]}$', fontsize=fs)
        self.ax[0].set_ylabel(r'Brain pdf', fontsize=fs)
        self.ax[0].set_xlim(self._run.A_lim[0],self._run.A_lim[1])
        self.ax[0].tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax[0].tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax[0].tick_params('both', length=12, width=2., which='major',
                             right=True, top=False)
        self.ax[0].tick_params('both', length=6, width=2., which='minor',
                             right=True, top=False) 
        self.ax[0].xaxis.set_minor_locator(MultipleLocator(0.2))
        self.ax[0].xaxis.set_major_locator(MultipleLocator(1.))
        self.ax[0].yaxis.set_minor_locator(MultipleLocator(0.1))
        self.ax[0].yaxis.set_major_locator(MultipleLocator(0.5))  
        self.ax[0].axes.get_yaxis().set_visible(False)

        self.ax[1].set_xlabel(r'$\rm{log}\ \tau\ \rm{[s]}$', fontsize=fs)
        self.ax[1].set_xlim(self._run.logtau_lim[0],self._run.logtau_lim[1])
        self.ax[1].tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax[1].tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax[1].tick_params('both', length=12, width=2., which='major',
                             right=True, top=False)
        self.ax[1].tick_params('both', length=6, width=2., which='minor',
                             right=True, top=False) 
        self.ax[1].xaxis.set_minor_locator(MultipleLocator(0.2))
        self.ax[1].xaxis.set_major_locator(MultipleLocator(1.))
        self.ax[1].axes.get_yaxis().set_visible(False)       
       
    def retrieve_data(self):
        fpath = ('./../OUTPUT_FILES/RUNS/' + self._run.subdir
                 + 'most_likely_pars.csv')
        self.df_mcmc = pd.read_csv(fpath, header=0, low_memory=False)

        #self.tau, self.tau_min, self.tau_max, self.A, self.B, self.C, self.chi2 = np.loadtxt(
        #  fpath, skiprows=1, delimiter=',', usecols=(1,2,3,4,7,10,16), unpack=True)

        fpath = ('./../OUTPUT_FILES/RUNS/' + self._run.subdir
                 + 'estimated_A_tau_B.csv')
        self.df_est = pd.read_csv(fpath, header=0, low_memory=False)
        
        
        #self.tau_est, self.A_est, self.B_est = np.loadtxt(
        #  fpath, skiprows=1, delimiter=',', usecols=(1,2,4), unpack=True)          
        
        #self.tau_unc = np.divide(np.maximum(self.tau_min,self.tau_max),self.tau) * 100.
        
        #Load values by Julien.

    def trim_data(self):
        #Remove voxels which might have bad fits. i.e. arbitrarily high tau.
        #cond = ((self.tau <= 200.) & (self.A >= -0.8) & (self.A <= 0.8) & (self.chi2 <= 4.e6))
        #cond = ((self.A <= 0.5))
        #cond = ((self.chi2 <= 1.e6))
        cond = ((self.tau_unc <= 25.))
        self.tau_est, self.A_est, self.B_est =\
          self.tau_est[cond], self.A_est[cond], self.B_est[cond]        
        self.tau, self.A, self.B, self.C, self.chi2, self.tau_unc =\
          self.tau[cond], self.A[cond], self.B[cond], self.C[cond],\
          self.chi2[cond], self.tau_unc[cond]
        print len(self.tau_est), len(self.tau), len(self.tau_unc)

    def plot_models(self):
        N_bins = 100. #Same number of steps in tau guesses when estimating pars.

        add_curves(
          self.ax[0],self.df_est['A_est'].values,(self._run.A_lim[0],
          self._run.A_lim[1]),N_bins,c[0])
        add_curves(
          self.ax[0],self.df_mcmc['A'].values,(self._run.A_lim[0],
          self._run.A_lim[1]),N_bins,c[1])
        add_curves(
          self.ax[0],self.df_mcmc['A'].values - self.df_mcmc['A_l'].values,
          (self._run.A_lim[0],self._run.A_lim[1]),N_bins,c[2])

        add_curves(
          self.ax[1],np.log10(self.df_est['tau_est'].values),
          (self._run.logtau_lim[0],self._run.logtau_lim[1]),N_bins,c[0])
        add_curves(
          self.ax[1],np.log10(self.df_mcmc['tau'].values),
          (self._run.logtau_lim[0],self._run.logtau_lim[1]),N_bins,c[1])
        add_curves(
          self.ax[1],np.log10(self.df_mcmc['tau'].values - self.df_mcmc['tau_l'].values),
          (self._run.logtau_lim[0],self._run.logtau_lim[1]),N_bins,c[2])
        
    def make_legend(self):
        self.ax[0].plot(
          [np.nan], [np.nan], marker='None', ls='-', lw=12., color=c[0],
          label=r'$\chi^2$')
        self.ax[0].plot(
          [np.nan], [np.nan], marker='None', ls='-', lw=12., color=c[2],
          label=r'MCMC: 16%')
        self.ax[0].plot(
          [np.nan], [np.nan], marker='None', ls='-', lw=12., color=c[1],
          label=r'MCMC: 50%')
        self.ax[0].legend(
          frameon=False, fontsize=fs, numpoints=1, labelspacing=0.2, loc=1)           

    def manage_output(self):
        if self._run.save_fig:
            fpath = os.path.join(
              './../OUTPUT_FILES/RUNS/' + self._run.subdir, 'FIGURES/Fig_hist.png')
            plt.savefig(fpath, format='png')
        if self._run.show_fig:
            plt.show()
        plt.close()
            
    def make_plot(self):
        self.set_fig_frame()
        self.retrieve_data()
        #self.trim_data()
        self.plot_models()
        self.make_legend()
        self.manage_output()

