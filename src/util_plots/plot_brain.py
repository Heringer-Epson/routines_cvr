#!/usr/bin/env python
import os, sys
sys.path.append(os.path.join(os.path.dirname(__file__), '..', 'lib'))

import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import nibabel as nib
import matplotlib.colors as colors

var2col_est = {'A': (2,)}
var2col_proc = {'A': (4,), 'A_unc': (4,5,6)}

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
        self.make_plot()
        
    def retrieve_data(self):
        
        #Load raw data.
        fpath = os.path.join('./../data_test/normal1/', 'CVR_raw_scaled.nii')
        img = nib.load(fpath)
        
        #Load brain mask -- mask out the voxels that fall outside the brain.
        fpath = os.path.join('./../data_test/normal1/', 'CVRmask.nii')
        mask_obj = nib.load(fpath)  
        mask = mask_obj.get_data().astype(bool)
        self.mask = mask     

        if self.mode == 'est':
            fpath = ('./../OUTPUT_FILES/RUNS/' + self._run.subdir
                     + 'estimated_A_tau_B.csv')
            cols = var2col_est[self.var]
            vals = np.loadtxt(
              fpath, skiprows=1, delimiter=',', usecols=cols, unpack=True)

        elif self.mode == 'proc':
            fpath = ('./../OUTPUT_FILES/RUNS/' + self._run.subdir
                     + 'most_likely_pars.csv')
            cols = var2col_proc[self.var]
            vals = np.loadtxt(
              fpath, skiprows=1, delimiter=',', usecols=cols, unpack=True)
        if vals.shape[0] == 3:
            qtty = np.abs(vals[0])
            unc = np.abs(np.maximum(vals[1],vals[2]))
            vals = np.divide(unc,qtty) * 100. #(in percentage)

        print 'Median of ploted values:', np.median(vals)
        print 'Maximum of ploted values:', max(vals)
        
        #Initialize parameter arrays and reshape the data.
        self.vals = np.zeros(mask.shape, dtype=float)        
        self.vals[mask] = vals
    
    def treat_data(self):
        if self.var == 'A':
            self.cmap = mpl.cm.get_cmap('PuOr')
            self.lims, self.linthresh, self.linscale = (-2.,2.), 0.05, 0.7
            self.norm = colors.SymLogNorm(
              self.linthresh,self.linscale,vmin=self.lims[0],vmax=self.lims[1])

        if self.var == 'A_unc':
            #self.lims = (0.,10.)
            self.lims = (1.,100.)
            self.cmap = mpl.cm.get_cmap('YlOrRd')
            #self.norm = colors.Normalize(vmin=self.lims[0], vmax=self.lims[1])
            self.norm = colors.LogNorm(vmin=self.lims[0], vmax=self.lims[1])

        #Trim data array.
        self.vals[self.vals <= self.lims[0]] = self.lims[0]
        self.vals[self.vals >= self.lims[1]] = self.lims[1]

    def loop_draw(self):

        N_cols = int(np.sqrt(self.vals.shape[-1])) + 1
        N_rows = self.vals.shape[-1] // N_cols + 1
    
        #Set frame.
        self.fig, ax = plt.subplots(N_rows, N_cols)
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
                if idx < self.vals.shape[-1] - 2:
                    data_2D = self.vals[:, :, idx]

                    #im = ax[i][j].imshow(
                      #data_2D.transpose(), interpolation='nearest', 
                    #  data_2D.transpose(), interpolation='sinc', 
                    #  norm=colors.SymLogNorm(
                    #  self.linthresh, self.linscale,
                    #  vmin=self.lims[0], vmax=self.lims[1]),
                    #  cmap=cmap, origin='upper', aspect='auto')            

                    im = ax[i][j].imshow(
                      data_2D.transpose(), interpolation='sinc', 
                      norm=self.norm,cmap=self.cmap, origin='upper', aspect='auto') 

                    out = np.ma.masked_array(self.mask[:, :, idx],self.mask[:, :, idx] == True)
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
