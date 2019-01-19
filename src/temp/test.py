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

class Test_Case(object):
    """
    Description:
    ------------
    Test case. Learning how to read .nii data.
    
    Parameters:
    -----------
    show_fig : ~bool
        Whether the created figure is to be displayed or not.
    save_fig : ~bool
        Whether the created figure is to be saved or not.
     
    Outputs:
    --------
    ./../../OUTPUT_FILES/FIGURES/test_mask.pdf
    
    References:
    -----------
    https://www.youtube.com/watch?v=9ffUQo2mF6w
    """   
        
    def __init__(self, show_fig, save_fig):
        
        self.show_fig = show_fig
        self.save_fig = save_fig

        fig = plt.figure(figsize=(16,8))
        self.ax1 = fig.add_subplot(131)
        self.ax2 = fig.add_subplot(132)
        self.ax3 = fig.add_subplot(133)
        
        self.img = None
        self.mask = None
        self.data = None
        self.mask_data = None
        self.unmasked_data = None
        self.run_test()

    def set_fig_frame(self):        
        x_label = r'coord $x_1$'
        y_label = r'coord $x_2$'
        self.ax1.set_xlabel(x_label, fontsize=fs)
        self.ax1.set_ylabel(y_label, fontsize=fs)
        self.ax2.set_xlabel(x_label, fontsize=fs)
        self.ax2.set_ylabel(y_label, fontsize=fs)
        
        self.ax1.tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax1.tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax1.minorticks_off()
        self.ax1.tick_params(
          'both', length=8, width=1., which='major', direction='in')
        self.ax1.tick_params(
          'both', length=4, width=1., which='minor', direction='in')
        self.ax1.xaxis.set_ticks_position('both')
        self.ax1.yaxis.set_ticks_position('both')             
        self.ax1.get_xaxis().set_visible(False)
        self.ax1.get_yaxis().set_visible(False)
        
        self.ax2.tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax2.tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax2.minorticks_off()
        self.ax2.tick_params(
          'both', length=8, width=1., which='major', direction='in')
        self.ax2.tick_params(
          'both', length=4, width=1., which='minor', direction='in')
        self.ax2.xaxis.set_ticks_position('both')
        self.ax2.yaxis.set_ticks_position('both')   
        self.ax2.get_xaxis().set_visible(False)
        self.ax2.get_yaxis().set_visible(False)

        self.ax3.tick_params(axis='y', which='major', labelsize=fs, pad=8)      
        self.ax3.tick_params(axis='x', which='major', labelsize=fs, pad=8)
        self.ax3.minorticks_off()
        self.ax3.tick_params(
          'both', length=8, width=1., which='major', direction='in')
        self.ax3.tick_params(
          'both', length=4, width=1., which='minor', direction='in')
        self.ax3.xaxis.set_ticks_position('both')
        self.ax3.yaxis.set_ticks_position('both')   
        self.ax3.get_xaxis().set_visible(False)
        self.ax3.get_yaxis().set_visible(False)

    def load_data(self):
        data_fpath = os.path.join('./../../data_test/normal1/', 'CVR_raw_scaled.nii')
        mask_fpath = os.path.join('./../../data_test/normal1/', 'CVRmask.nii')
        self.img = nib.load(data_fpath)
        self.mask = nib.load(mask_fpath)
    
    def read_test_data(self):
        
        #from nibabel.testing import data_path
        #fpath = os.path.join(os.path.join(data_path, 'example4d.nii.gz'))

        #print self.img.shape
        #print self.img.header
        #print dir(self.img)
        #print dir(self.img.dataobj)
        #print self.img.dataobj
        #vol1 = self.img.dataobj[10, 10, 10, 10]
        #print dir(vol1)
        #print header['cal_max']

        self.data = self.img.get_data()
        self.mask_data = self.mask.get_data().astype(bool)

    def compute_std(self):
        #Calculate std deviation and re-assign it to a shape that include
        #outside the brain.
        std_data_masked = np.std(self.data[self.mask_data], axis=1)
        #print std_data_masked.shape
        
        self.unmasked_data = np.zeros(self.mask_data.shape, dtype=std_data_masked.dtype)
        self.unmasked_data[self.mask_data] = std_data_masked
        #print self.unmasked_data.shape        

    def plot_mask(self):
        self.ax1.imshow(
          self.data[:, :, 30, 90], interpolation='nearest', cmap='viridis')
        self.ax1.set_title('EPI Image', fontsize=fs)
        self.ax2.imshow(
          1. - self.mask_data[:, :, 30], interpolation='nearest', cmap='Reds')
        self.ax2.set_title('Mask Image', fontsize=fs)        
        self.ax3.imshow(
          self.unmasked_data[:, :, 30], interpolation='nearest', cmap='bone')
        self.ax3.set_title('std', fontsize=fs)  

    def reshape_data(self):

        #This will reshape all of the original data in a 2D array,
        #where the first entry is the combinate x,y,z dimensions and the
        #second in the time.
        data_2D = self.data.reshape(
          np.prod(self.data.shape[:-1]), self.data.shape[-1])
        #print data_2D.shape

        #More effectively, one may use the mask to not only ignore the
        #voxels outside the brain, but also reshape the data at the same time.
        #print self.data[self.mask_data].shape
        
        print self.data.shape
        print self.data[30,30,30,10]
        
        
    def manage_output(self):
        if self.save_fig:
            fpath = os.path.join('./../OUTPUT_FILES/FIGURES/', 'test_mask.pdf')
            plt.savefig(fpath, format='pdf')
        if self.show_fig:
            plt.show()         
        
    def run_test(self):
        self.set_fig_frame()
        self.load_data()
        self.read_test_data()
        self.compute_std()
        self.plot_mask()
        self.reshape_data()
        self.manage_output()
        
if __name__ == '__main__':  
    Test_Case(show_fig=False, save_fig=True)


