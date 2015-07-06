import numpy as np
from pylab import *

"""
    SCRIPT NAME:
    radxfer.py

    AUTHOR:
    Greg Blumberg
    wblumberg@ou.edu

    DESCRIPTION:
    The primary argument for this code is the opd, which are the optical depths 
    set at the wave number grid (it's a 2D array of dimensions (layers, wnums)).

    This code also takes a wavenumber grid, a temperature profile, a surface temperature,
    a surface emissivity, and a Boolean designating whether or not the calculation is
    upwelling or downwelling and performs the radiative transfer calculations.

    The temperature profile should be of the dimension (layers) and can come from the
    tape7_reader.py file.

    This code was based off of Dave Turner's (NSSL) IDL code.
"""

def planck(wnum, temp):
    c1 = 1.1910427e-5 #mW/m2.sr.cm-1
    c2 = 1.4387752 #K/cm
    r = (c1 * np.power(wnum,3)) / ( np.exp( c2*(wnum/temp)) - 1.)
    return r

def rt(wnum, temp, opd, sfc_t=None, sfc_e=None, upwelling=False):
    opd = np.asarray(opd, dtype=np.float64)
    wnum = np.asarray(wnum, dtype=np.float64)
    temp = np.asarray(temp, dtype=np.float64)
    if upwelling is False:
        sfce = 0.
        k_start = len(opd)-1
        k_end = -1
        k_step = -1
    else:
        k_start = 0
        k_end = len(opd)
        k_step = 1
    rad = np.zeros(len(wnum))

    for k in np.arange(k_start, k_end, k_step):
        trans = np.asarray(np.exp(-1. * opd[k,:]), dtype=np.float64)
        if upwelling is True:
            b_boundary = planck(wnum, temp[k+1])
        else:
            b_boundary = planck(wnum, temp[k])
        b_avg = planck(wnum, (temp[k] + temp[k+1])/2.)
        rad = rad * trans + (1.-trans) * (b_boundary+2.*(b_avg-b_boundary) * (1./opd[k,:] - trans / (1.-trans)))
    return rad

if __name__ == '__main__':
    print rt([2,3,4],[5,6,7],[[2,3,4],[1,2,3]])
    
