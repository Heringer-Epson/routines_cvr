import sys, os
import cPickle
import numpy as np
import lib.stats as stats
from scipy.interpolate import interp1d
from lib import stats
from lib import data_handling 

def fv(_var):
    return ',' + str(format(_var, '.6f'))

def vars2line1(_idx, _vars):
    _A, A_l, A_u, _tau, tau_l, tau_u = _vars
    return fv(_A) + fv(A_l) + fv(A_u) + fv(_tau) + fv(tau_l) + fv(tau_u)

def vars2line2(_idx, _vars):
    _A, _tau, _B = _vars
    return fv(_A) + fv(_tau) + fv(_B)

class Compute_Likelihoods(object):
    """
    Code Description
    ----------    
    Loads a .pkl file containing models for each cell of the parameter space
    and assign (compute) a likelihood to each model. Likelihoods are computed
    using a Bayesian method. The likelihood depends on the combined probability
    of the points belonging to a region (in time) of the datasample. This
    region is defined as an input parameter (self._run.region_to_fit). For a
    given model, the probability of each point is computed by assuming that the
    data follows a (normal) gaussian, centered at that data point and whose
    sigma is determined by the noise. The code marginalises the probability
    over the 'baseline' variable.
    
    Parameters
    ----------
    _run : object
      instance of input_params containing the data in each case of study.

    Outputs:
    --------
    ./../OUTPUT_FILES/RUNS/' + self._run.subdir, 'most_likely_A_tau.csv     
    ./../OUTPUT_FILES/RUNS/' + self._run.subdir, 'most_likely_A_tau_B.csv     
    """
    def __init__(self, _run):
        self._run = _run        
        self.models = data_handling.get_convolved_models(self._run)        

    def iterate_voxels(self):
        #Here, the suffixes 1 and 2 denote marginalised and non-marginalised,
        #respectively.
        
        L = {}
        L['parspace'] = self._run.parspace
        L['parspace_reduced'] = self._run.parspace_reduced
        columns = zip(*self._run.parspace_reduced)
        As, taus = np.array(columns[0]), np.array(columns[1])
        cond = data_handling.region2cond(self._run)
        
        fpath1 = os.path.join(
          './../OUTPUT_FILES/RUNS/' + self._run.subdir, 'most_likely_A_tau.csv')
        fpath2 = os.path.join(
          './../OUTPUT_FILES/RUNS/' + self._run.subdir, 'most_likely_A_tau_B.csv')
        with open(fpath1, 'w') as out1, open(fpath2, 'w') as out2:
            out1.write('voxel,A,A_unc_l,A_unc_u,tau,tau_unc_l,tau_unc_u')
            out2.write('voxel,A,tau,B')
         
            #for idx_space in range(self._run.signal_s.shape[0]):
            #    bold_voxel = self._run.signal_s[idx_space,:]

            for idx_space in [3523]: #fewer voxels. 15821
                bold_voxel = self._run.signal_ns[idx_space,:]
                y_obs = bold_voxel[cond]
                y_obs_unc = self._run.pCO2_noise[cond]   
                
                ln_L1, ln_L2 = [], []
                for (A,tau) in self._run.parspace_reduced:
                    
                    #[0:] is so that the first data point is not include. By
                    #definition, its convolution is always zero.
                    ln_L = np.array([stats.compute_L(
                      y_obs,self.models[
                        data_handling.pars2label(A,tau,B)][cond][0:],y_obs_unc[0:])
                        for B in self._run.B])
                    ln_L1.append(stats.marginalize(ln_L, self._run.B))
                    ln_L2.extend(ln_L)
                    
                ln_L1, ln_L2 = np.array(ln_L1), np.array(ln_L2)
                L['i_' + str(idx_space)] = ln_L1
                L['i_red_' + str(idx_space)] = ln_L2
                
                outvars = stats.get_contour_uncertainties(As,taus,ln_L1)
                print outvars
                
                #outvars = stats.plot_contour(
                #  None,As,taus,ln_L1,'b',add_max=True,show_fig=False)
                #line = vars2line1(idx_space,outvars)
                #out1.write('\n' + str(idx_space) + line)
                #print outvars

                idx_max = ln_L2.argmax()
                A_max, tau_max, B_max = self._run.parspace[idx_max]
                line = vars2line2(idx_space,(A_max,tau_max,B_max))
                out2.write('\n' + str(idx_space) + line)

        fpath = os.path.join(
          './../OUTPUT_FILES/RUNS/' + self._run.subdir, 'likelihoods.pkl')
        with open(fpath, 'w') as output:
            cPickle.dump(L, output, cPickle.HIGHEST_PROTOCOL)    

    def run_task(self):
        self.iterate_voxels()
