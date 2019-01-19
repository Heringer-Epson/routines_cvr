import sys, os, time
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import cPickle
import numpy as np
from scipy.interpolate import interp1d
import core_funcs as cf
import stats
import data_handling 

def fv(_var):
    return ',' + str(format(_var, '.6f'))
def vars2line1(_idx, _vars):
    _A, A_l, A_u, _tau, tau_l, tau_u = _vars
    return fv(_A) + fv(A_l) + fv(A_u) + fv(_tau) + fv(tau_l) + fv(tau_u)

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
    """
    #@profile
    def __init__(self, _run):
        print 'Calculating likelihoods...'
        L = {}

        A_qE = (_run.A_step / 2.)**2.
        
        fpath_S = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'PICKLES/smooth.pkl'
        fpath_M = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'PICKLES/models.pkl'
        fpath_out = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'most_likely_A_tau.csv'
        N_trials = 0
        with open(fpath_S, 'r') as fS, open(fpath_M, 'r') as fM,\
          open(fpath_out, 'w') as out:
            S, M = cPickle.load(fS), cPickle.load(fM)
            out.write('voxel,A,A_unc_l,A_unc_u,tau,tau_unc_l,tau_unc_u')

            columns = zip(*_run.parspace_Atau)
            As, taus = np.array(columns[0]), np.array(columns[1])
            cond = data_handling.region2cond(_run,S['time'])
            N_voxels = S['signal_ns'].shape[0]
            ln_B_range = np.log(_run.B[-1] - _run.B[0])
      
            start_time = time.time()
            print '    N voxels = ' + str(N_voxels)
            for idx_space in range(N_voxels):
            #for idx_space in [3523]: #fewer voxels. 15821
                print idx_space
                N_trials += 1
                bold_voxel = S['signal_ns'][idx_space,:]
                signal_unc = S['signal_noise'][idx_space,:]
                pCO2_unc = S['pCO2_noise']
                
                ln_L = []
                for (A,tau) in _run.parspace_Atau:
                    label = data_handling.pars2label(A,tau,0.)
                    
                    unc = np.sqrt(signal_unc**2. + (A * pCO2_unc)**2.) #improve.                  

                    ln_L_of_B = np.array([
                      cf.compute_L(bold_voxel[cond],M[
                      data_handling.pars2label(A,tau,B)][cond],unc[cond])
                      for B in _run.B])                     
                    
                    ln_L_of_B_clean, max_base = stats.treat_array(ln_L_of_B)
                    ln_L_margin = (cf.marg(ln_L_of_B_clean, _run.B) + max_base - ln_B_range)                    
                                        
                    ln_L.append(ln_L_margin)
   
                ln_L = np.array(ln_L)
                L['likelihood_list_' + str(idx_space)] = ln_L
                                
                outvars = stats.get_contour_uncertainties(As,taus,ln_L,A_qE)
                line = vars2line1(idx_space,outvars)
                
                out.write('\n' + str(idx_space) + line)
                delta_time = time.time() - start_time

            data_handling.save_pickle(_run.subdir, 'likelihoods.pkl', L)
                
        print '    Run took ', format(delta_time, '.1f'), 's'
        #print '    Approx ', format(delta_time / float(N_voxels), '.3f'), ' s/voxel'
        print '    Approx ', format(delta_time / float(N_trials), '.3f'), ' s/voxel'
        print '    Done.\n'
