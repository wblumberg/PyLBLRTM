import numpy as np

"""
    SCRIPT NAME:
    apodizer.py

    SCRIPT AUTHOR:
    Greg Blumberg

    DATE MADE: 
    9/25/2013

    SUMMARY:
    This script utilizies the output from the convolve2aeri.py script
    in order to apodize the spectra (i.e. remove the higher frequency
    portions of the spectra...with the optical depth values,
    this removes ODs that are less than 0.)

"""

def apodizer(spectrum):
    n = len(spectrum)
    imd = n / 2
    apod = apodize_norton_beer(n, imd)
    new_spectrum = np.fft.ifft(np.fft.fft(spectrum)*apod)
    return new_spectrum

#n - length of the apodization function (double sized ifg)
#md - index to point where opd = MOPD (for single sized ifg)
def apodize_norton_beer(n, idx_md):
    beer = np.zeros(n)
    beer[0] = 1.
    for i in np.arange(1,n/2):
        if (i <= idx_md):
            beer[i] = (1-((i-1)/float(idx_md))**2)**2
        else:
            beer[i] = 0.

    if (n % 2 == 0):
        beer[n/2:n] = beer[0:n/2][::-1]
    else:
        beer[n/2:n] = beer[0:(n/2)+1][::-1]
    return beer
    
