import sys, os
import cPickle
import numpy as np
import lib.stats as stats
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
from lib.make_A_tau_model import Make_Model
from util_plots.plot_voxel import Plot_Signal
from util_plots.plot_example_contour import Plot_Contour
from lib import stats

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
        
        self.models = {}

    def get_convolved_models(self):
        fpath = os.path.join(
          './../OUTPUT_FILES/RUNS/' + self._run.subdir, 'models.pkl')        
        if os.path.isfile(fpath): 
            with open(fpath, 'r') as f:
                self.models = cPickle.load(f)
        else:
            raise ValueError(
              'File ' + fpath + ' has not been created yet. Please '
              + 'run master.py with flag run_models=True.')

        #Check if the parameter space targeted is the same as the one used to
        #compute likelihoods.
        if not np.array_equal(self.models['parspace'],self._run.parspace):
            raise ValueError(
              'Frror: parspace defined in input_params.py does not match the '
              + 'array saved in models.pkl. Please re-run master.py with flag '
              + 'run_models=True')             

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
                cond = ((self._run.time >= 0.) & (self._run.time <= 1.e6))                    
         
            #for idx_space in range(self._run.signal_s.shape[0]):
            #    bold_voxel = self._run.signal_s[idx_space,:]

            for idx_space in range(1): #fewer voxels.
                bold_voxel = self._run.signal_ns[15821,:] #3523
                #bold_voxel = self._run.signal_ns[3523,:] #

                ln_L_coll = []
                for p, (A,tau,B) in enumerate(self._run.parspace):
                    label = 'conv_' + str(p)
                    
                    y_obs = bold_voxel[cond]
                    y_obs_unc = self._run.pCO2_noise[cond]            
                    
                    ln_L = stats.compute_L(y_obs,self.models[label][cond],y_obs_unc)
                    #ln_L_of_A_tau.append(ln_L_test)                    
                    #print _ln_L, ln_L_test
                
                #print self.models['conv_10'], self._run.parspace[10]
                #Plot_Signal(self._run.time, bold_voxel, self.models['conv_5000'], self._run.pCO2)
                Plot_Signal(self._run.time, bold_voxel, self.models['conv_5000'] + B_fit, self._run.pCO2)
                
                #Plot_Signal(self._run.time, bold_voxel, self.models[best_label], self._run.pCO2)
                
                ln_L_of_A_tau = np.array(ln_L_of_A_tau)

           
                outvars = stats.plot_contour(
                  self.ax,As,taus,ln_L_of_A_tau,'b',add_max=True,show_fig=False)
                line = vars2line(idx_space,outvars)
                out.write(line)
                         
                idx_max = ln_L_of_A_tau.argmax()
                best_label = 'conv_' + str(p)
                print self._run.parspace[idx_max]
                
                #print self._run.B[idx_space]
                
                Plot_Signal(self._run.time, bold_voxel, self.models[best_label] + B_fit, self._run.pCO2)
                #Plot_Contour(self.models['parspace'], ln_L_of_A_tau)
                                   
                #Check zero at the beginning of all models.
                #model sigma to include pCO2 and signal noise.
                #Make plotting routine to show one example of cvr best fit and contours.
                #make histogram of best fits.
                #Models are missing a baseline value....
                #Use pkl run for everything.


    '''
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
                cond = ((self._run.time >= 0.) & (self._run.time <= 1.e6))                    
         
            #for idx_space in range(self._run.signal_s.shape[0]):
            #    bold_voxel = self._run.signal_s[idx_space,:]

            for idx_space in range(1): #fewer voxels.
                bold_voxel = self._run.signal_ns[15821,:] #3523
                #bold_voxel = self._run.signal_ns[3523,:] #

                B_fit = self._run.B[idx_space]
                B_unc = self._run.B_unc[idx_space]
                marg_range = np.linspace(
                  B_fit - 2. * B_unc, B_fit + 2. * B_unc + 1.e-4, self._run.B_trials)

                ln_L_of_A_tau = []
                for p, (A,tau) in enumerate(self._run.parspace):
                    label = 'conv_' + str(p)
                    
                    y_obs = bold_voxel[cond]
                    y_obs_unc = self._run.pCO2_noise[cond]

                    #ln_L_marg = np.array([stats.compute_L(
                    #  y_obs,self.models[label][cond] + _B,y_obs_unc) for _B in
                    #  marg_range])
                    #_ln_L = stats.marginalize(ln_L_marg, marg_range, B_fit, B_unc)
                    #ln_L_of_A_tau.append(_ln_L)                    
                    
                    #ln_L_test = stats.compute_L(y_obs,self.models[label][cond] + B_fit,y_obs_unc)
                    #ln_L_of_A_tau.append(ln_L_test)                    
                    #print _ln_L, ln_L_test
                
                print self._run.parspace[5000]
                #print self.models['conv_10'], self._run.parspace[10]
                #Plot_Signal(self._run.time, bold_voxel, self.models['conv_5000'], self._run.pCO2)
                Plot_Signal(self._run.time, bold_voxel, self.models['conv_5000'] + B_fit, self._run.pCO2)
                
                #Plot_Signal(self._run.time, bold_voxel, self.models[best_label], self._run.pCO2)
                
                ln_L_of_A_tau = np.array(ln_L_of_A_tau)

           
                outvars = stats.plot_contour(
                  self.ax,As,taus,ln_L_of_A_tau,'b',add_max=True,show_fig=False)
                line = vars2line(idx_space,outvars)
                out.write(line)
                         
                idx_max = ln_L_of_A_tau.argmax()
                best_label = 'conv_' + str(p)
                print self._run.parspace[idx_max]
                
                #print self._run.B[idx_space]
                
                Plot_Signal(self._run.time, bold_voxel, self.models[best_label] + B_fit, self._run.pCO2)
                #Plot_Contour(self.models['parspace'], ln_L_of_A_tau)
                                   
                #Check zero at the beginning of all models.
                #model sigma to include pCO2 and signal noise.
                #Make plotting routine to show one example of cvr best fit and contours.
                #make histogram of best fits.
                #Models are missing a baseline value....
                #Use pkl run for everything.
    '''
    def run_task(self):
        self.get_convolved_models()
        self.iterate_voxels()

