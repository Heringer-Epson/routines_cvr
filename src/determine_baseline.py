import numpy as np

class Fit_Baseline(object):
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
        self.B, self.B_unc = [], []

    def run_fit(self):
        
        cond = ((self._run.time >= self._run.rest[0])
                & (self._run.time <= self._run.rest[1]))
        time_rest = self._run.time[cond]
        
        for idx_space in range(self._run.signal_ns.shape[0]):
            bold_voxel = self._run.signal_ns[idx_space,:]
            out_B, out_B_unc = np.polyfit(time_rest,bold_voxel[cond],0,cov=True)
            self.B.append(out_B[0])
            self.B_unc.append(np.sqrt(1. / out_B_unc[0][0]))
             
    def run_task(self):
        self.run_fit()
        return np.array(self.B), np.array(self.B_unc)
               
