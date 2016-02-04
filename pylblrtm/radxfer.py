import numpy as np
from pylab import *
import sys

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
    r = (c1 * np.power(wnum,3)) / ( np.exp( c2*wnum/temp) - 1.)
    return r

def rt(wnum, temp, opd, zenith_angle=0, sfc_t=None, sfc_e=None, upwelling=False, debug=False):
    opd = np.asarray(opd, dtype=np.float64) * (1./np.cos(np.radians(zenith_angle)))
    wnum = np.asarray(wnum, dtype=np.float64)
    temp = np.asarray(temp, dtype=np.float64)
    if upwelling is False:
        sfce = 0.
        k_start = len(opd)-1
        k_end = -1
        k_step = -1
    else:
        sfce = sfc_e
        if sfce < 0 or sfce > 1:
            print "Error: the surface emissivity is outside of [0,1]"
            sys.exit()        
        if sfc_t is None:
            # if no surface temperature is given, use the temperature of the lowest level
            sfct = temp[0]
        else:
            sfct = sfc_t
        sfc_b = planck(wnum, sfct)
        k_start = 0
        k_end = len(opd)
        k_step = 1

    rad = np.zeros(len(wnum))

    # use the linear in tau approach
    if sfce > 0:
        rad = sfc_b * sfce + rt(wnum, temp, opd) * (1.-sfce)
    # B_close is the Planck function for the temperature at the edge of the layer 
    # B_far | --> | B_close
    # 
    for k in np.arange(k_start, k_end, k_step)[::-1]:
        trans = np.asarray(np.exp(-1. * opd[k,:]), dtype=np.float64) # Compute the transmissivity of the layer
        if upwelling is True:
            # Compute the b_boundary from the bottom layer up
            b_close = planck(wnum, temp[k])
            layer_to_inst = np.exp(-1. * np.sum(opd[:k,:], axis=0))
        else:
            # Compute the b_boundary from the top of the layer down
            b_close = planck(wnum, temp[k+1])
            layer_to_inst = np.exp(-1. * np.sum(opd[k+1:,:], axis=0))
        b_avg = (planck(wnum, temp[k]) + planck(wnum, temp[k+1]))/2. # b_close and b_far
        b_eff = b_close + 2*(b_avg - b_close)*((1./opd[k,:]) - (trans/(1.-trans)))
        if debug == True:
            print "Temperature of top and bottom layer:", temp[k], temp[k+1]
            print "Planck top and bottom layer:", planck(wnum, temp[k]), planck(wnum, temp[k+1])
            print "b_avg:", b_avg
            print "Temperature of top and bottom layer:", temp[k], temp[k+1]
            print "b_close:", b_close 
            print "b_avg:", b_avg
            print "b_eff:", b_eff
            print "Optical Depth of Current Layer:", opd[k,:]
            print "Terms of the RT for this layer:", (1-np.exp(-1.*opd[k,:])), b_avg, layer_to_inst
            print "Calculation:", (1-np.exp(-1.*opd[k,:]))*b_eff*layer_to_inst 
        rad = rad*trans + (1.-trans) * b_eff
    return rad

if __name__ == '__main__': 
    #plot(np.arange(100,3000,1), planck(np.arange(100,3000,1), 300))
    #show()
    # problem 1
    # down=40.69 and up=57.16
    ods = [[0.002, 0.008, 0.032, 0.125, 0.5]]
    #ods = [[0.2, 0.8, 3.2, 12.5, 50.0]]
    ods = np.asarray(ods).T
    temp = [230, 244, 258, 272, 286, 300]
    emiss = 1
    wnum = [1200]
    #wnum = [800]

    print "TEMP:", temp
    print "Optical Depths:", ods
    print "WNUM:", wnum
    upw_rad = rt(wnum, temp, ods, zenith_angle=60, sfc_t=310, sfc_e=0.9, upwelling=True)
    print "Upwelling Radiance:", upw_rad
    print
    dwn_rad = rt(wnum, temp, ods, zenith_angle=60) 
    print "Downwelling Radiance:", dwn_rad

  



 
