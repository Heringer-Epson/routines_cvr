import numpy as np
from scipy.signal import savgol_filter

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
        self._run = _run
        self.signal_n, self.signal_ns, self.signal_sn_unc = None, None, None
        self._pCO2_s, self.pCO2_noise = None, None

    def normalize_data(self):
        #Compute the mean signal in time for each voxel and divide each voxel's
        #signal by its mean in time. Do the latter for all time steps.
        self.signal_n = np.copy(self._run.signal)
        bold_median = np.median(self._run.signal, axis=1)
        for idx_time in range(self._run.signal.shape[1]):
            self.signal_n[:,idx_time] = np.divide(
              self._run.signal[:,idx_time],bold_median) * 100.

    def smooth(self):
        self.signal_ns = savgol_filter(
          self.signal_n,self._run.smoothing_window,3)
        self._pCO2_s = savgol_filter(
          self._run.pCO2,self._run.smoothing_window,3)        
    
    def estimate_noise(self):
        self.pCO2_noise = compute_rms(
          self._run.pCO2,self._pCO2_s,self._run.time,10.)
        self.signal_noise = compute_rms(
          self.signal_n,self.signal_ns,self._run.time,10.)
    
    def run_task(self):
        self.normalize_data()
        self.smooth()
        self.estimate_noise()
        return self.signal_n, self.signal_ns, self.signal_sn_unc,\
               self._pCO2_s, self.pCO2_noise
               
