import sys, os, time
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import cPickle
import numpy as np
import core_funcs as cf
import stats
import emcee
import corner 
import matplotlib.pyplot as plt
from multiprocessing import Pool
from functools import partial

#http://dfm.io/emcee/current/user/line/

#@profile
def lnprob(theta, y, yerr, t, ts, I, tau_l, tau_u, A_l, A_u, B_l, B_u, C_l,
           C_u, D_l, D_u):
    tau, A, B, C, D = theta
    if (
      tau_l < tau < tau_u and A_l < A < A_u and B_l < B < B_u
      and C_l < C < C_u and D_l < D < D_u): #Priors.
        x_conv = I(t, tau)
        return cf.lnlike(A, B, C, D, ts, x_conv, y, yerr)
    else:
        return -np.inf  

def pars2line(idx,tau,A,B,C,D,chi2):
    def S(p):
        if np.isnan(p):
            return ',NaN'
        else:
            return ',' + str(format(p, '.4f'))   
    return (
      '\n' + str(idx) + S(tau[0]) + S(tau[1]) + S(tau[2]) + S(A[0]) + S(A[1])
      + S(A[2]) + S(B[0]) + S(B[1]) + S(B[2]) + S(C[0]) + S(C[1]) + S(C[2]) 
      + S(D[0]) + S(D[1]) + S(C[2]) + S(chi2))

def apply_MCMC(S, I, thetas, t, ts, tau_lim, A_lim, B_lim, C_lim, D_lim, idx_space):
    theta_est = thetas[idx_space]
    if not np.isnan(theta_est[2]):
        y = S['signal_ns'][idx_space,:]
        yerr = S['signal_noise'][idx_space,:]
        ndim, nwalkers = 5, 100
        pos = [theta_est + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
        sampler = emcee.EnsembleSampler(
          nwalkers, ndim, lnprob, args=(y, yerr, t, ts, I, tau_lim[0], 
          tau_lim[1], A_lim[0], A_lim[1], B_lim[0], B_lim[1], C_lim[0],
          C_lim[1], D_lim[0], D_lim[1]))
        sampler.run_mcmc(pos, 500)
        samples = sampler.chain[:, 50:, :].reshape((-1, ndim))
        tau_mcmc, A_mcmc, B_mcmc, C_mcmc , D_mcmc = map(
          lambda v: (v[1], v[2]-v[1], v[1]-v[0]),
          zip(*np.percentile(samples, [16, 50, 84], axis=0)))
        #Compute chi2
        x_conv = I(t, tau_mcmc[1])
        chi2 = cf.chi2(A_mcmc[1], B_mcmc[1], C_mcmc[1], D_mcmc[1], ts, x_conv, y, yerr)
        
    #fig = corner.corner(
    #  samples, labels=["$\\tau$", "$A$", "B"],
    #  truths=[theta_est[0], theta_est[1], theta_est[2]])
    #plt.show()

    else:
        val = (np.nan,np.nan,np.nan,np.nan,np.nan)
        tau_mcmc, A_mcmc, B_mcmc, C_mcmc, chi2 = val, val, val, val    
    return pars2line(idx_space, tau_mcmc, A_mcmc, B_mcmc, C_mcmc, D_mcmc, chi2)

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
          fpath_est, skiprows=1, delimiter=',', usecols=(1,2,4), unpack=True)

        #Load interpolator for computing pCO2 convolutions.
        fpath_I = './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'PICKLES/tau_interp.pkl'
        fI = open(fpath_I, 'r')
        self.I = cPickle.load(fI)
        fI.close()
        
    def compute_best_pars(self):

        print 'Calculating likelihoods...'       
        fpath_out = './../OUTPUT_FILES/RUNS/' +self. _run.subdir + 'most_likely_pars.csv'

        with open(fpath_out, 'w') as out:
            out.write(
              'voxel,tau,tau_l,tau_u,A,A_l,A_u,B,B_l,B_u,C,C_l,C_u,D,D_l,'
              + 'D_u,chi2')

            N_voxels = self.S['signal_ns'].shape[0]      
            t = self.S['time']
            ts = t / self._run.t_pivot
            C_est = np.zeros(len(self.tau_est))
            D_est = np.zeros(len(self.tau_est))

            thetas = np.transpose(
              np.array([self.tau_est, self.A_est, self.B_est, C_est, D_est]))
            
            start_time = time.time()

            #Parallel.
            output = []
            MCMC_given_idx = partial(
              apply_MCMC, self.S, self.I, thetas, t, ts, self._run.tau_lim,
              self._run.A_lim, self._run.B_lim, self._run.C_lim, self._run.D_lim)
            pool = Pool(5)
            output += pool.map(MCMC_given_idx,range(N_voxels))
            #output += pool.map(MCMC_given_idx,range(2000))
            pool.close()
            pool.join()
            for line in output:
                out.write(line) 

            #Non-parallel.
            #for idx_space in range(10):
            #    out.write(apply_MCMC(
            #      self.S, self.I, thetas, t, ts, self._run.tau_lim,
            #       self._run.A_lim, self._run.B_lim, self._run.C_lim,
            #       self._run.D_lim, idx_space))
        print '    Run took ', format(time.time() - start_time, '.1f'), 's'
        print '    Done.\n'

    def run_model(self):
        self.retrieve_data()
        self.compute_best_pars()
