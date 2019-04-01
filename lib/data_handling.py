import sys, os
import cPickle
import numpy as np

def save_pickle(_subdir, _fname, _D):
    fpath = './../OUTPUT_FILES/RUNS/' + _subdir + 'PICKLES/' + _fname
    if os.path.isfile(fpath): 
        os.remove(fpath)
    with open(fpath, 'w') as output:
        cPickle.dump(_D, output, cPickle.HIGHEST_PROTOCOL)             

def pars2label(_A,_tau,_B):
    p1 = str(format(_A, '.6f'))
    p2 = str(format(_tau, '.6f'))
    p3 = str(format(_B, '.6f'))
    return 'm_' + p1 + '_' + p2 + '_' + p3
    
def region2cond(_run,_t):
    if _run.region_to_fit == 'rest':
        cond = ((_t >= _run.rest[0]) & (_t <= _run.rest[1]))
    elif _run.region_to_fit == 'step':
        cond = ((_t >= _run.step[0]) & (_t <= _run.step[1]))
    elif _run.region_to_fit == 'ramp':
        cond = ((_t >= _run.ramp[0]) & (_t <= _run.ramp[1]))
    elif _run.region_to_fit == 'norest':
        cond = ((_t >= _run.norest[0]) & (_t <= _run.norest[1]))
    elif _run.region_to_fit == 'all':
        cond = ((_t >= -1.e4) & (_t <= 1.e4))
    return cond


    
