import sys, os, time
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import cPickle
import numpy as np
import core_funcs as cf
import stats
import emcee
import corner 
import matplotlib.pyplot as plt

#http://dfm.io/emcee/current/user/line/

#@profile
def lnprob(theta, y, yerr, t, I):
    tau, A, B = theta
    if -10. < A < 10. and 10. < B < 200. and 1. < tau < 100.: #Priors.
        x_conv = I(t, tau)
        return cf.lnlike(A, B, x_conv, y, yerr)
    else:
        return -np.inf    

def pars2line(idx,tau,A,B):
    def S(p):
        if np.isnan(p):
            return ',NaN'
        else:
            return ',' + str(format(p, ',f'))   
    return (
      '\n' + str(idx) + S(tau[0]) + S(tau[1]) + S(tau[2]) + S(A[0]) + S(A[1])
      + S(A[2]) + S(B[0]) + S(B[1]) + S(B[2]))
      
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
    def __init__(self, _run):
        self._run = _run
        self.tau_est, self.A_est, self.B_est = None, None, None
        self.run_model()

    def retrieve_data(self):

        #Load observational data.
        fpath_S = './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'PICKLES/smooth.pkl'
        fS = open(fpath_S, 'r')
        self.S = cPickle.load(fS)
        fS.close()        
        
        #Load par estimative.
        fpath_est = './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'estimated_A_tau_B.csv'
        self.tau_est, self.A_est, self.B_est = np.loadtxt(
          fpath_est, skiprows=1, delimiter=',', usecols=(1,2,3), unpack=True)

        #Load interpolator for computing pCO2 convolutions.
        fpath_I = './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'PICKLES/tau_interp.pkl'
        fI = open(fpath_I, 'r')
        self.I = cPickle.load(fI)
        fI.close()
        
    def compute_best_pars(self):

        print 'Calculating likelihoods...'
        fpath_out = './../OUTPUT_FILES/RUNS/' +self. _run.subdir + 'most_likely_pars.csv'

        with open(fpath_out, 'w') as out:
            out.write('voxel,tau,tau_l,tau_u,A,A_l,A_u,B,B_l,B_u')

            N_voxels = self.S['signal_ns'].shape[0]      
            t = self.S['time']
            
            start_time = time.time()
            for idx_space in range(N_voxels):
            #for idx_space in range(10):
                print idx_space

                #theta_est = np.array(
                #  [self.tau_est[idx_space], self.A_est[idx_space],
                #  self.B_est[idx_space]]) 

                theta_est = np.array(
                  [self.tau_est, self.A_est,
                  self.B_est]) 
                
                if not np.isnan(theta_est[2]):

                    y = self.S['signal_ns'][idx_space,:]
                    yerr = self.S['signal_noise'][idx_space,:]
                    

                    
                    start_time = time.time()
                    
                    ndim, nwalkers = 3, 100
                    pos = [theta_est + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
                    
                    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob,
                                                    args=(y, yerr, t, self.I))
                    sampler.run_mcmc(pos, 500)
                    
                    samples = sampler.chain[:, 50:, :].reshape((-1, ndim))

                    tau_mcmc, A_mcmc, B_mcmc = map(
                      lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
                      zip(*np.percentile(samples, [16, 50, 84], axis=0)))

                    #print tau_mcmc, A_mcmc, B_mcmc
                    #print theta_est
                    
                    #fig = corner.corner(
                    #  samples, labels=["$\\tau$", "$A$", "B"],
                    #  truths=[theta_est[0], theta_est[1], theta_est[2]])
                    #plt.show()
                else:
                    val = (np.nan,np.nan,np.nan)
                    tau_mcmc, A_mcmc, B_mcmc = val, val, val
                    
                out.write(pars2line(idx_space, tau_mcmc, A_mcmc, B_mcmc))

                                                    
        delta_time = time.time() - start_time
    
        print '    Run took ', format(delta_time, '.1f'), 's'
        #print '    Approx ', format(delta_time / float(N_voxels), '.6f'), ' s/voxel'
        print '    Done.\n'

    def run_model(self):
        self.retrieve_data()
        self.compute_best_pars()
