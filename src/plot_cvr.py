import os
import numpy as np
import nibabel as nib
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib.ticker import MultipleLocator

mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['mathtext.fontset'] = 'stix'
mpl.rcParams['font.family'] = 'STIXGeneral'
fs = 24.

class Plot_Cvr(object):
    """
    Description:
    ------------
    Test case. Plot cvr as a function of time.
    
    Parameters:
    -----------
    show_fig : ~bool
        Whether the created figure is to be displayed or not.
    save_fig : ~bool
        Whether the created figure is to be saved or not.
     
    Outputs:
    --------
    ./../../OUTPUT_FILES/FIGURES/X.pdf
    
    References:
    -----------
    https://www.youtube.com/watch?v=9ffUQo2mF6w
    """   
        
    def __init__(self, show_fig, save_fig):
        
        self.show_fig = show_fig
        self.save_fig = save_fig

        fig = plt.figure(figsize=(10,10))
        self.ax = fig.add_subplot(111)
        
        self.img = None
        self.data = None
        self.pco2 = None
        self.bold_median = None
        self.make_plot()

    def set_fig_frame(self):        
        x_label = r'coord $x_1$'
        y_label = r'coord $x_2$'
        self.ax.set_xlabel(x_label, fontsize=fs)
        self.ax.set_ylabel(y_label, fontsize=fs)
        self.ax.set_xlabel(x_label, fontsize=fs)
        self.ax.set_ylabel(y_label, fontsize=fs)
        
        self.ax.tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax.tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax.minorticks_off()
        self.ax.tick_params(
          'both', length=8, width=1., which='major', direction='in')
        self.ax.tick_params(
          'both', length=4, width=1., which='minor', direction='in')
        self.ax.xaxis.set_ticks_position('both')
        self.ax.yaxis.set_ticks_position('both')             
        self.ax.get_xaxis().set_visible(False)
        self.ax.get_yaxis().set_visible(False)
        
    def load_data(self):
        data_fpath = os.path.join('./../data_test/normal1/', 'CVR_raw_scaled.nii')
        #data_fpath = os.path.join('./../data_test/normal1/', 'CVR.nii')
        mask_fpath = os.path.join('./../data_test/normal1/', 'CVRmask.nii')
        pco2_fpath = os.path.join('./../data_test/normal1/', 'pco2_0.1D')
        
        self.img = nib.load(data_fpath)
        mask_obj = nib.load(mask_fpath)
        self.pco2 = np.loadtxt(pco2_fpath, usecols=(0,), unpack=True)
        
        mask = mask_obj.get_data().astype(bool)
        self.data = self.img.get_data()
        self.data = self.data[mask].astype(float) #Data only includes voxels within the brain.
        #print self.data.shape
        #print pco2.shape

    def process_bold(self):
        '''When using the mean, the procedure below leads to values that agree
        to 0.01 with the values obtained by Julien and stored in 'CVR_raw_scaled.nii'. 

        Q: why the mean in time? Doesn't seem reliable nor intuitive. I'd
        think that the mean in volume gives a better representation of how the
        whole brain is respoding to signals, making it easier to identify the
        signal from sick voxels.
        '''
                  
        #Compute the mean signal in time for each voxel.
        self.bold_median = np.median(self.data, axis=1)
        #At each time, fivide each voxel's signal by its mean in time.
        print self.bold_median.shape
        for i in range(self.data.shape[1]):
            self.data[:,i] /= self.bold_median

        
        #print self.data.shape
        #print self.data[1,:]
        #self.bold_median = np.median(self.data, axis=0)
        #print self.bold_median
        #self.data = np.divide(self.data,self.bold_median)
        #print self.data.shape

        #print self.bold_median.shape, self.pco2.shape
        #print self.data[0,:].shape
        #print self.data[100,:]

    def plot_data(self):
        t = np.linspace(0, 100, len(self.pco2))
        self.ax.plot(t, self.data[100], ls='-', color='k')
        self.ax.plot(t, self.pco2, ls='-', color='r')

    def manage_output(self):
        if self.save_fig:
            fpath = os.path.join('./../OUTPUT_FILES/FIGURES/', 'X.pdf')
            plt.savefig(fpath, format='pdf')
        if self.show_fig:
            plt.show()         
        
    def make_plot(self):
        self.set_fig_frame()
        self.load_data()
        self.process_bold()
        self.plot_data()
        self.manage_output()
        
if __name__ == '__main__':  
    Plot_Cvr(show_fig=True, save_fig=False)


