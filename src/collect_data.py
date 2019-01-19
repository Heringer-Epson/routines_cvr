import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import cPickle
import numpy as np
import data_handling 

class Collect_Data(object):
    """
    Description:
    ------------
    Read data from .nii file which is passed as an input parameters. 
    
    Parameters:
    -----------
    _run : object
      instance of input_params containing the data in each case of study.

    Outputs:
    --------
    Creates './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'master.pkl'
    """   

    def __init__(self, _run):
        import nibabel as nib
        M = {}
        
        #Load the data.
        data_fpath = os.path.join(_run.dirpath, 'CVR.nii')
        mask_fpath = os.path.join(_run.dirpath, 'CVRmask.nii')
        pco2_fpath = os.path.join(_run.dirpath, 'pco2_0.1D')
        
        self.img = nib.load(data_fpath)
        mask_obj = nib.load(mask_fpath)
        M['pCO2'] = np.loadtxt(pco2_fpath, usecols=(0,), unpack=True)
        
        mask = mask_obj.get_data().astype(bool)
        signal = self.img.get_data()
        M['signal'] = signal[mask].astype(float) #Data only includes voxels within the brain.

        #make a time array.
        N_time = M['signal'].shape[1]
        dt = _run.time_step
        duration = N_time * dt
        M['time'] = np.arange(dt, duration + 1.e-3, dt)      

        data_handling.save_pickle(_run.subdir, 'data.pkl', M)

