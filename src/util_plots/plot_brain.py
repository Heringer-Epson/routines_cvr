#!/usr/bin/env python
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import matplotlib as mpl
import nibabel as nib
import matplotlib.colors as colors
from numpy.ma import masked_array
from nilearn.plotting.cm import _cmap_d


var2col_est = {'tau': (1,), 'A': (2,), 'A_unc': (3,), 'B': (3,), 'B_unc': (4,)}
var2col_proc = {'tau': (1,), 'tau_unc': (1,2,3), 'A': (4,), 'A_unc': (4,5,6),
                'B': (7,), 'A_unc': (7,8,9)}

#https://matplotlib.org/users/colormapnorms.html
class Plot_Brain(object):
    """
    Description:
    ------------
    TBW.
    
    Parameters
    ----------
    _run : object
      instance of input_params containing the data in each case of study.
    
    Outputs:
    --------
    TBW.
    """        
    def __init__(self, _run, var, mode):
        print 'Plotting brain maps...'
                
        self._run = _run
        self.var, self.mode = var, mode
        
        self.lims, self.linthresh, self.linscale = None, None, None
        self.img, self.mask, self.data, self.vals = None, None, None, None
        self.fig, self.norm, self.cmap = None, None, None
        self.anat = None
        self.make_plot()
        
    def retrieve_data(self):
        
        #Load raw data.
        #fpath = os.path.join(self._run.dirpath, 'CVR_raw_scaled.nii')
        #img = nib.load(fpath)
        
        #Load brain mask -- mask out the voxels that fall outside the brain.
        fpath = os.path.join(self._run.dirpath, 'CVRmask.nii')
        mask_obj = nib.load(fpath)  
        mask = mask_obj.get_data().astype(bool)
        self.mask = mask

        #Initialize variables.
        self.A = np.zeros(mask.shape, dtype=float)
        self.tau = np.zeros(mask.shape, dtype=float)

        #Load relevant file. Estimated or processed values.
        if self.mode == 'est':
            fpath = ('./../OUTPUT_FILES/RUNS/' + self._run.subdir
                     + 'estimated_A_tau_B.csv')
            self.df = pd.read_csv(fpath, header=0, low_memory=False)
            self.A[mask] = self.df['A_est'].values
            self.tau[mask] = self.df['tau_est'].values
        elif self.mode == 'proc':
            fpath = ('./../OUTPUT_FILES/RUNS/' + self._run.subdir
                     + 'most_likely_pars.csv')
            self.df = pd.read_csv(fpath, header=0, low_memory=False)
            self.A[mask] = self.df['A'].values
            self.tau[mask] = self.df['tau'].values

        '''
        #Tests load data from .nii output file.
        proc_fpath = './../OUTPUT_FILES/RUNS/' + self._run.subdir + 'proc.nii'
        proc_obj = nib.load(proc_fpath)
        proc_data = proc_obj.get_data()
        self.tau = proc_data[:,:,:,1]
        self.A = proc_data[:,:,:,4]
        '''
  
    def treat_data(self):
        if self.var == 'A':
            self.cmap = _cmap_d['cold_white_hot']
            self.lims = (-.7,.7)
            self.norm = colors.Normalize(vmin=self.lims[0],vmax=self.lims[1])
            self.A[self.A <= self.lims[0]] = self.lims[0]
            self.A[self.A >= self.lims[1]] = self.lims[1]

        if self.var == 'tau':
            self.cmap = mpl.cm.get_cmap('plasma_r')
            self.tau = np.log10(self.tau)
            self.lims = (0.,2.)
            self.norm = colors.Normalize(vmin=self.lims[0],vmax=self.lims[1])
            self.tau[self.tau <= self.lims[0]] = self.lims[0]
            self.tau[self.tau >= self.lims[1]] = self.lims[1]

    def loop_draw(self):

        N_cols = int(np.sqrt(self.A.shape[-1])) + 1
        N_rows = self.A.shape[-1] // N_cols + 1
    
        #Set frame.
        self.fig, ax = plt.subplots(N_rows, N_cols, figsize=(15,12))
        self.fig.patch.set_facecolor('k')
        plt.subplots_adjust(wspace=0, hspace=0, left=0, right=1., top=1, bottom=0)  

        #Make drawings.
        self.cmap.set_bad(color='black')
        for i in range(N_rows):
            for j in range(N_cols):
                idx = i * N_cols + j
                ax[i][j].set_facecolor('k')
                ax[i][j].set_frame_on(False)
                ax[i][j].axes.get_xaxis().set_visible(False)
                ax[i][j].axes.get_yaxis().set_visible(False)
                if idx < self.A.shape[-1] - 2:              
                    A_2D = self.A[:, :, idx].transpose()      
                    tau_2D = self.tau[:, :, idx].transpose()        

                    if self.var == 'A':
                        im = ax[i][j].imshow(
                          A_2D, interpolation='spline36',norm=self.norm, 
                          cmap=self.cmap, origin='upper', aspect='auto') 

                    elif self.var == 'tau':


                        im = ax[i][j].imshow(
                          tau_2D, interpolation='spline36', 
                          norm=self.norm,cmap=self.cmap,
                          origin='upper', aspect='auto') 

                        #tau_2D_Apos = masked_array(tau_2D,A_2D<0.) #Masks out!
                        #tau_2D_Aneg = masked_array(tau_2D,A_2D>=0.) #Masks out!                        #im = ax[i][j].imshow(
                        #  tau_2D_Apos, interpolation='nearest', 
                        #  norm=self.norm,cmap=mpl.cm.get_cmap('summer_r'),
                        #  origin='upper', aspect='auto') 
                        #im = ax[i][j].imshow(
                        #  tau_2D_Aneg, interpolation='nearest', 
                        #  norm=self.norm,cmap=mpl.cm.get_cmap('cool'),
                        #  origin='upper', aspect='auto') 

                    out = np.ma.masked_array(
                      self.mask[:, :, idx],self.mask[:, :, idx] == True)
                    im = ax[i][j].imshow(
                      out.transpose(), interpolation='nearest', 
                      cmap='Greys_r', origin='upper', aspect='auto')  

    def manage_output(self):
        if self._run.save_fig:
            fname = 'Fig_map_' + self.var + '_' + self.mode + '.pdf'
            fpath = os.path.join(
              './../OUTPUT_FILES/RUNS/' + self._run.subdir, 'FIGURES/' + fname)
            plt.savefig(fpath, format='pdf', facecolor=self.fig.get_facecolor())
        if self._run.show_fig:
            plt.show()
        plt.close()

    def make_plot(self):
        self.retrieve_data()
        self.treat_data()
        self.loop_draw()
        self.manage_output()
