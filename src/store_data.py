import sys, os
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import numpy as np
import pandas as pd

class Store_Data(object):
    """
    Description:
    ------------
    TBW.
    
    Parameters:
    -----------
    _run : object
      instance of input_params containing the data in each case of study.

    Instructions:
    -------------
    https://www.youtube.com/watch?v=9ffUQo2mF6w

    Outputs:
    --------
    Creates './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'master.pkl'
    """   

    def __init__(self, _run):
        import nibabel as nib

        mask_fpath = os.path.join(_run.dirpath, 'CVRmask.nii')
        mask_obj = nib.load(mask_fpath)
        
        #Collected object properties from mask. Header and affine will be
        #copied into the new .nii file.
        mask_data = mask_obj.get_data().astype(bool)
        mask_header = mask_obj.get_header()
        mask_affine = mask_obj.get_affine()
        unmasked_canvas = np.zeros(mask_data.shape, dtype=float)

        #Get processed data.
        fpath = './../OUTPUT_FILES/RUNS/' + _run.subdir + '/most_likely_pars.csv'
        df_proc = pd.read_csv(fpath, header=0, low_memory=False)
        tau = df_proc['tau'].values        
        proc = df_proc.values
        N_cols = proc.shape[1]
        
        #Create canvas for 4D matrix (volume + processed data)
        unmasked_proc_data = np.zeros((
          mask_data.shape[0],mask_data.shape[1],mask_data.shape[2],N_cols))       
        
        #There probably is a more efficient way to do this.
        for i in range(N_cols):
            unmasked_proc_data[:,:,:,i][mask_data] = proc[:,i]
        
        mask_header.set_data_dtype(unmasked_proc_data.dtype)
        out_img = nib.Nifti1Image(unmasked_proc_data, mask_affine, mask_header)
        out_img.to_filename('./../OUTPUT_FILES/RUNS/' + _run.subdir + '/proc.nii')

