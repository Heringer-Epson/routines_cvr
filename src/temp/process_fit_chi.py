import sys, os
import numpy as np
import nibabel as nib
from scipy.optimize import curve_fit
from scipy.integrate import simps, quad
from scipy.interpolate import interp1d
from util_plots.plot_voxel import Plot_Signal

def make_exp(tau, t_i, t_f):
    #norm = 1. / (tau * (np.exp(-t_i / tau) - np.exp(-t_f / tau)))
    norm = 1.
    def _exp_func(t):
        return norm * np.exp(-t / tau)
    return _exp_func        

def compute_convolution(_pCO2_func,_time):
    def model(t,A,tau):
        exp_func = make_exp(tau, _time[0], _time[-1])
        cond = (_time <= t)        
        _t = _time[cond]
        a = _pCO2_func(_t)
        b = exp_func(- (t - _t) / tau)
        _expected_signal = simps(np.multiply(a,b), _t)
        return _expected_signal
    return model

def scaling_model(_pCO2):
    def model(t,A):
        expected_signal = A * _pCO2
        return expected_signal
    return model    

class Process_Chi(object):
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
        self.time_fine = None

    def normalize_data(self):
        #Compute the mean signal in time for each voxel and divide each voxel's
        #signal by its mean in time. Do the latter for all time steps.
        bold_median = np.median(self._run.signal, axis=1)
        bold_median_s = np.median(self._run.signal_s, axis=1)
        for idx_time in range(self._run.signal.shape[1]):
            self._run.signal[:,idx_time] /= (bold_median / 100.)
            self._run.signal_s[:,idx_time] /= (bold_median / 100.)
        
    def iterate_voxels(self):
        
        if self.mode == 'step':
            cond = ((self._run.time >= self._run.step[0])
                    & (self._run.time <= self._run.step[1]))
            self.time_fine = np.arange(
              self._run.time[cond][0], self._run.time[cond][-1], 0.01)
        
        #for idx_space in range(self._run.data.shape[0]):
        for idx_space in range(1): #fewer voxels.
            #bold_voxel = self._run.signal[idx_space,:]
            bold_voxel = self._run.signal_s[15821,:]
            #bold_voxel = self._run.signal_s[3523,:]
        
            pCO2_func = interp1d(
              self._run.pCO2[cond],self._run.time[cond],bounds_error=False,
              fill_value=(self._run.pCO2[cond][0],self._run.pCO2[cond][1]))
                        
            _model = np.vectorize(compute_convolution(pCO2_func,self.time_fine))
                        
            popt, pcov = curve_fit(_model, self._run.time[cond], bold_voxel[cond],p0=[2.5,.5])
            print popt
            #signal_fit = [_model(_t,popt[0],popt[1]) for _t in self._run.time]           
            #signal_fit = [_model(_t,2.5,1.) for _t in self._run.time]           
            
            
            #NO CONV.
            #_model = scaling_model(self._run.pCO2)
            #popt, pcov = curve_fit(_model, self._run.time, bold_voxel,p0=[1.])
            
            #signal_fit = _model(self._run.time,popt[0])           
            
            
            #if self._run.show_fig:
            #    Plot_Signal(self._run.time, bold_voxel, signal_fit, self._run.pCO2)


    def run_task(self):
        self.normalize_data()
        self.iterate_voxels()

            
if __name__ == '__main__':  
    pass
    #Process_Chi()
