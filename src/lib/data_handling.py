import sys, os
import cPickle
import numpy as np
import nibabel as nib

class Read_Data(object):
    """
    Description:
    ------------
    Read data from .nii file which is passed as an input parameters. 
    
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
        self._signal, self._pCO2 = None, None

    def load_data(self):
        data_fpath = os.path.join(self._run.dirpath, 'CVR.nii')
        mask_fpath = os.path.join(self._run.dirpath, 'CVRmask.nii')
        pco2_fpath = os.path.join(self._run.dirpath, 'pco2_0.1D')
        
        self.img = nib.load(data_fpath)
        mask_obj = nib.load(mask_fpath)
        self._pCO2 = np.loadtxt(pco2_fpath, usecols=(0,), unpack=True)
        
        mask = mask_obj.get_data().astype(bool)
        self._signal = self.img.get_data()
        self._signal = self._signal[mask].astype(float) #Data only includes voxels within the brain.

    def make_time_array(self):
        N_time = self._signal.shape[1]
        dt = self._run.time_step
        duration = N_time * dt
        self._time_array = np.arange(dt, duration + 1.e-3, dt)      
        
    def run_task(self):
        self.load_data()
        self.make_time_array()
        return self._signal, self._pCO2, self._time_array


def get_convolved_models(_run):
    fpath = os.path.join('./../OUTPUT_FILES/RUNS/' + _run.subdir, 'models.pkl')        
    if os.path.isfile(fpath): 
        with open(fpath, 'r') as f:
            models = cPickle.load(f)
    else:
        raise ValueError(
          'File ' + fpath + ' has not been created yet. Please run master.py '\
          + 'with flag run_models=True.')

    #Check if the parameter space targeted is the same as the one used to
    #compute likelihoods.
    if not np.array_equal(models['parspace'],_run.parspace):
        raise ValueError(
          'Frror: parspace defined in input_params.py does not match the '
          + 'array saved in models.pkl. Please re-run master.py with flag '
          + 'run_models=True')    

    return models

def pars2label(_A,_tau,_B):
    p1 = str(format(_A, '.6f'))
    p2 = str(format(_tau, '.6f'))
    p3 = str(format(_B, '.6f'))
    return 'm_' + p1 + '_' + p2 + '_' + p3
    
    
    
    
    
    
