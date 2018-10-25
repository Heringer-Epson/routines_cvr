#!/usr/bin/env python
import numpy as np
import pandas as pd
from scipy.integrate import simps

"""
Set of functions which can be used to calculate likelihoods.
"""       

def clean_array(inp_array):
    """Limit the number of magnitudes (to 20 orders below maximum) in the input
    array and then normalize it. This is helpful for dealing with likelihoods
    whose log spans ~1000s of mags. Input should be an array of the natural
    log of the values in interest. output array is in linear scale.
    """
    inp_array[inp_array < max(inp_array) -20.] = max(inp_array) - 20.     
    inp_array = inp_array - min(inp_array) 
    inp_array = np.exp(inp_array)
    inp_array = inp_array / sum(inp_array)
    return inp_array 

def get_contour_levels(inp_array, contour):
    """Given an input array that is normalized (i.e. cumulative histogram adds
    to one), return the values in that array that correspond to the confidence
    limits passed in 'contour'.
    """
    _L_hist = sorted(inp_array, reverse=True)
    _L_hist_cum = np.cumsum(_L_hist)

    _L_hist_diff = [abs(value - contour) for value in _L_hist_cum]
    diff_at_contour, idx_at_contour = min((val,idx) for (idx,val)
                                          in enumerate(_L_hist_diff))
    #Check if contour placement is too coarse (>10% level).
    if diff_at_contour > 0.1:
        UserWarning(str(contour * 100.) + '% contour not constrained.')	
    return _L_hist[idx_at_contour]

def plot_contour(ax, x, y, z, c, add_max=True, show_fig=False):

    x_best = x[np.argmax(z)]
    y_best = y[np.argmax(z)]

    try:
        contour_list = [0.95, 0.68, 0.] 
        z = clean_array(z)
        _x, _y = np.unique(x), np.unique(y)       
        X = x.reshape(len(_x),len(_y))
        Y = y.reshape(len(_x),len(_y))
        qtty = z.reshape((len(_x), len(_y)))
        levels = [get_contour_levels(z, contour) for contour in contour_list]

        ax.contourf(X, Y, qtty, levels[0:2], colors=c, alpha=0.4)	 
        cs = ax.contourf(X, Y, qtty, levels[1:3], colors=c, alpha=0.6)	 

        if add_max:
            ax.plot(x_best, y_best, ls='None', marker='+',color=c, markersize=10.)
        ax.plot([np.nan], [np.nan], color=c, ls='-', lw=15., marker='None')

        #Estimate parameter uncertainty from ellipsis path.
        p = cs.collections[0].get_paths()[0]
        v = p.vertices
        x_cont, y_cont = v[:,0], v[:,1]    
        x_unc = (max(x_cont) - x_best, x_best - min(x_cont))
        y_unc = (max(y_cont) - y_best, y_best - min(y_cont))
    except:
        x_unc, y_unc = (np.nan,np.nan), (np.nan,np.nan)
    
    return x_best, x_unc[0], x_unc[1], y_best, y_unc[0], y_unc[1]

def compute_L(y_obs, y_pred, noise):
    """..."""
    ln_L = np.log(1. / (2 * np.pi * noise**2.)) - (y_obs - y_pred)**2. / (2. * noise**2.) 
    return np.sum(ln_L)

def marginalize(ln_L_of_theta, theta):
    
    ln_L_of_theta[ln_L_of_theta < max(ln_L_of_theta) -20.] = max(ln_L_of_theta) - 20.
    
    norm = simps(np.ones(len(theta)), theta)
    
    max_ln_L = max(ln_L_of_theta)
    ln_l_aux = ln_L_of_theta - max_ln_L
    ln_out = np.log(simps(np.exp(ln_l_aux),theta)) + max_ln_L - np.log(norm)    
    return ln_out
    
    
    
    
    


