import sys, os
import cPickle
import numpy as np
import lib.stats as stats
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from util_plots.plot_voxel import Plot_Signal
from util_plots.plot_example_contour import Plot_Contour
from lib import stats
from lib.data_handling import get_convolved_models
from lib.data_handling import pars2label

def norm_dist(x, mu, sigma):
    norm = 1. / np.sqrt(2. * np.pi * sigma**2.)
    exp_factor = np.exp(-(x - mu)**2. / (2. * sigma**2.))
    return norm * exp_factor
    
def fv(_var):
    return str(format(_var, '.4f')) + ','

def vars2line(_idx, _vars):
    _A, A_l, A_u, _tau, tau_l, tau_u = _vars
    _line = ('\n' + fv(_A) + fv(A_l) + fv(A_u) + fv(_tau) + fv(tau_l) + fv(tau_u))
    return _line[0:-1]

class Compute_Likelihoods(object):
    """
    Code Description
    ----------    
    TBW.

    Parameters
    ----------
    TBW.
    """
    def __init__(self, _run, mode):
        self._run = _run
        self.mode = mode
        
        self.fig = plt.figure(figsize=(10,10))
        self.ax = self.fig.add_subplot(111)
        
        self.models = get_convolved_models(self._run)        

    def iterate_voxels(self):
        columns = zip(*self._run.parspace)
        As, taus = np.array(columns[0]), np.array(columns[1])
        
        fpath = os.path.join(
          './../OUTPUT_FILES/RUNS/' + self._run.subdir, 'most_likely_A_tau.csv')
        with open(fpath, 'w') as out:
            out.write('voxel,A,A_unc_l,A_unc_u,tau,tau_unc_l,tau_unc_u')
    
            if self.mode == 'step':
                #cond = ((self._run.time >= self._run.step[0])
                #        & (self._run.time <= self._run.step[1]))

                #cond = ((self._run.time >= self._run.ramp[0])
                #        & (self._run.time <= self._run.ramp[1]))

                cond = ((self._run.time >= 0.) & (self._run.time <= 1.e6))                    
         
            #for idx_space in range(self._run.signal_s.shape[0]):
            #    bold_voxel = self._run.signal_s[idx_space,:]

            for idx_space in range(1): #fewer voxels.
                                
                #bold_voxel = self._run.signal_ns[15821,:] #3523
                bold_voxel = self._run.signal_ns[3523,:] #
                y_obs = bold_voxel[cond]
                y_obs_unc = self._run.pCO2_noise[cond]   
                
                ln_L_coll = []
                for (A,tau) in self._run.parspace_reduced:
                    
                    ln_L_marg = np.array([stats.compute_L(
                      y_obs,self.models[pars2label(A,tau,B)][cond],y_obs_unc)
                      for B in self._run.B])
                    
                    if ((A == 0.25) & (tau == 1.)):
                        print '\n', ln_L_marg
                    
                    ln_L_coll.append(stats.marginalize(ln_L_marg, self._run.B))
                

                
                ln_L_coll = np.array(ln_L_coll)
                idx_max = ln_L_coll.argmax()
                A_max, tau_max = self._run.parspace_reduced[idx_max]
                best_label = pars2label(A_max, tau_max, 90.)                
                Plot_Signal(self._run.time, bold_voxel, self.models[best_label], self._run.pCO2)
                
           
                #outvars = stats.plot_contour(
                #  self.ax,As,taus,ln_L_of_A_tau,'b',add_max=True,show_fig=False)
                #line = vars2line(idx_space,outvars)
                #out.write(line)
                            
                #Check zero at the beginning of all models.
                #model sigma to include pCO2 and signal noise.
                #Make plotting routine to show one example of cvr best fit and contours.
                #make histogram of best fits.
                #Models are missing a baseline value....
                #Use pkl run for everything.
                #Marginalization gives significantly diff results. Check.
                #How to make plot for marginalized?

    #No marginalization.
    '''
    def iterate_voxels(self):
        columns = zip(*self._run.parspace)
        As, taus = np.array(columns[0]), np.array(columns[1])
        
        fpath = os.path.join(
          './../OUTPUT_FILES/RUNS/' + self._run.subdir, 'most_likely_A_tau.csv')
        with open(fpath, 'w') as out:
            out.write('voxel,A,A_unc_l,A_unc_u,tau,tau_unc_l,tau_unc_u')
    
            if self.mode == 'step':
                cond = ((self._run.time >= 0.) & (self._run.time <= 1.e6))                    

            for idx_space in range(1): #fewer voxels.
                bold_voxel = self._run.signal_ns[3523,:] #

                ln_L_coll = []
                for (A,tau,B) in self._run.parspace:
                    label = pars2label(A,tau,B)
                    y_obs = bold_voxel[cond]
                    y_obs_unc = self._run.pCO2_noise[cond]       
                    ln_L = stats.compute_L(y_obs,self.models[label][cond],y_obs_unc)
                    ln_L_coll.append(ln_L)
                
                ln_L_coll = np.array(ln_L_coll)
                idx_max = ln_L_coll.argmax()
                A_max, tau_max, B_max = self._run.parspace[idx_max]
                best_label = pars2label(A_max, tau_max, B_max)                
                Plot_Signal(self._run.time, bold_voxel, self.models[best_label], self._run.pCO2)
    '''
            


    def run_task(self):
        self.iterate_voxels()

