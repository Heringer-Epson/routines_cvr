import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import numpy as np
import cPickle
from scipy.signal import savgol_filter
import data_handling 

def compute_rms(y_data, y_smot, x_data, x_bin):
    """
    """                          
    def rms(_y_data, _y_smot):
        #Given a noisy and a smoothed data, compute an array of the
        #squared differences and take the square-root of its mean.
        #Used as a proxy of the noise.
        return np.sqrt(((_y_data - _y_smot)**2.).mean())

    #Compute the rms as a proxy of the noise of the flux point-wise.
    #Although point-wise, the noise of a given point is determined by
    #computing the rms including also the nearby points -- this prevents
    #funky values from being generated. In the loop below, for each point
    #'w' in the wavelength array, created a mini array containing the
    #nearby normalized and smoothed fluxes, which are then used as inputs
    #to the rms function.
    corr = 1.1
    y_rms = np.asarray([rms(
      y_data[(x_data >= x - x_bin) & (x_data <= x + x_bin)],
      y_smot[(x_data >= x - x_bin) & (x_data <= x + x_bin)])
      * corr for x in x_data])
    return y_rms

class Norm_Smooth(object):
    """
    Description:
    ------------
    TBW. 
    
    Parameters:
    -----------
    fpath : ~str
      Full path to the directory containing the MRI, pCO2 and mask data.

    Outputs:
    --------
    None
    """   
        
    def __init__(self, _run):
        print 'Smoothing the data...'
        S = {}
        fpath_D = './../OUTPUT_FILES/RUNS/' + _run.subdir + 'PICKLES/data.pkl'
        with open(fpath_D, 'r') as fD:
            D = cPickle.load(fD)
            
            #Copy a few relevant data from data.pkl into the smooth dict.
            S['time'], S['ts'], S['pCO2'] = D['time'], D['ts'], D['pCO2']
            try:
                S['grey_frac'], S['white_frac'], S['csf_frac'] =\
                  D['grey_frac'], D['white_frac'], D['csf_frac']
            except:
                pass
            
            #Compute the mean signal in time for each voxel and divide each voxel's
            #signal by its mean in time. Do the latter for all time steps.
            S['signal_n'] = np.zeros(shape=D['signal'].shape)
            bold_median = np.median(D['signal'], axis=1)
            for idx_time in range(D['signal'].shape[1]):
                S['signal_n'][:,idx_time] = np.divide(
                  D['signal'][:,idx_time],bold_median) * 100.

            #Smooth the signal.
            S['signal_ns'] = savgol_filter(S['signal_n'],_run.smoothing_window,3)
            S['pCO2_s'] = savgol_filter(S['pCO2'],_run.smoothing_window,3)        
        
            #Estimate noise
            S['pCO2_noise'] = compute_rms(S['pCO2'],S['pCO2_s'],S['time'],10.)
            S['signal_noise'] = np.zeros(
              shape=(D['signal'].shape[0],D['signal'].shape[1]))
            for idx_space in range(D['signal'].shape[0]):
                S['signal_noise'][idx_space,:] = compute_rms(
                S['signal_n'][idx_space,:],S['signal_ns'][idx_space,:],S['time'],10.)
        
            data_handling.save_pickle(_run.subdir, 'smooth.pkl', S)
        print '    Done.\n'

               
