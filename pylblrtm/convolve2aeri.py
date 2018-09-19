#import scipy.signal
import numpy as np
import sys
from . import panel_file

"""
    SCRIPT:
    convolve_to_aeri

    AUTHOR:
    Greg Blumberg
    wblumberg@ou.edu

    DESCRIPTION:
    This script takes an AERI radiance spectra and wavenumber grid from the
    LBLRTM monochromatic calculations and convolves it with the AERI 
    response function to produce an AERI radiance observation.

    It was designed after Dave Turner's IDL code.
"""

def convolve_to_aeri(wnum, spectra):

#ods = panel_file.panel_file('../idl_ex/out/ODdeflt_001', do_load_data=True)    

#wnum = ods.v
#spectra = ods.data1

#if True:   
    #Get minimum and maximum wavenumber values from the spectrum
    minv = np.min(wnum)
    maxv = np.max(wnum)

    #Ensure that the wavenumber (x) and the spectra (y) are 1D and are float64 arrays
    x = np.asarray(wnum, dtype=np.float64).flatten()
    y = np.asarray(spectra, dtype=np.float64).flatten()

    #Check for equal sizes in the wnum and spectra arrays
    if len(wnum) != len(spectra):
        print("ERROR: wnum and spectra are not the same number of elements.")
        sys.exit()

    #AERI wavenumber grid spacing
    AERI_dv = 15799./(2**15)
    AERI_OPD = 1./(2.*AERI_dv)

    #Step 1:
    #Apply tapers to the end of the radiance spectrum.  This will
    #result in the spectrum having 2^n + 1 points (which is optimum for FFT)
    #Find the mean wavenumber delta for this spectrum

    delx = np.mean(np.abs(np.diff(x)))
    
    #Step 2:
    #Taper the head (long wavelength side) of the spectrum
    npts = ((x[0] - 0)/delx) - 1 
    pts = ((x[0] - 0)/delx) - 1  #need number of points to the origin
    xx = np.arange(0,npts,1, dtype=np.float64) * delx 
    yy = np.zeros((npts+1))
    x = np.hstack((xx, x))
    y = np.hstack((yy, y))
    
    #Step 3:
    #Insert taper at the tail to get the proper numver of points in the spectrum (2^n + 1 points)
    for i in np.arange(0,101):
        if (2**i >= len(x)):
            break
    npts = 2**(i+1) - len(x) + 1
    xx = (np.arange(0, npts,1, dtype=np.float64) * delx) + delx + x[len(x) - 1]
    yy = np.zeros(len(xx))
    x = np.hstack((x, xx))
    y = np.hstack((y, yy))
    
    #Determine the size of the rolloff to apply; it should be no more than
    #100 cm-1, but it may need to be smaller?
    rolloffsize = np.min([np.max(x)-maxv, minv, 100.])
    tapersize = 20. #The amount of space in the spectrum to taper [cm-1]
    
    #Find the specral regions that require a roll-off
    v_rolloff1 = np.where((minv - rolloffsize <= x) & (x <= minv))[0]
    v_rolloff2 = np.where((maxv <= x) & (x <= maxv + rolloffsize))[0]
    
    len_v1 = len(v_rolloff1)
    len_v2 = len(v_rolloff2)
    #Apply the roll-off and then make it smooth for a small wavenumber region around the roll-off
    oldy = y #For debugging (in Dave's code)
    
    #Apply for the rolloff 1 points
    taper_function1 = (np.cos(np.arange(0,len_v1,1, dtype=np.float64)/(len_v1 - 1.) * np.pi - np.pi) + 1.)/2.
    idxs = np.where((minv <= x) & (x <= (minv+tapersize)))[0]
    len_idxs = len(idxs)
    y[v_rolloff1] = taper_function1 * np.mean(y[idxs])
    weights = np.arange(0, len_idxs,1, dtype=np.float32)/(len_idxs - 1.)
    y[idxs] = y[idxs] * weights + (1.-weights)*np.mean(y[idxs])
    
    #Apply for the rolloff 2 points
    taper_function2 = (np.cos(np.arange(0,len_v2,1, dtype=np.float64)/(len_v2 - 1.) * np.pi) + 1.)/2.
    idxs = np.where(((maxv - tapersize) <= x) & (x <= maxv))[0]
    len_idxs = len(idxs)
    y[v_rolloff2] = taper_function2 * np.mean(y[idxs])
    weights = np.arange(0, len_idxs,1, dtype=np.float32)/(len_idxs - 1.)
    y[idxs] = y[idxs] * weights + (1.-weights)*np.mean(y[idxs])
    
    #If the wavenumber res. is too coarse, then we need to zeropad the spectrum to allow us to interpolate
    #to multiple of the AERI res.
    if delx > 0.01:
        n = len(y)
        yfold = np.hstack((y, y[1:n-1][::-1])) #not n-2 because python doesn't include the next element in this syntax
        n = len(yfold)
        y_intemp = np.fft.fft(yfold)
        yyi = np.roll(np.asarray(y_intemp, dtype=np.float64), -1*n/2)

        #Want to zeropad the spectrum so we have 2^14 points 
        #Need figure out how many zeros to put at the ends of the interferogram
        #And we need to keep track of the factor that we are expanding the interferogram by so
        #we can mutiply the spectrum by it in a later step to put the energy back into it
        for i in np.arange(0,101):
            if 2**i >= len(x):
                break
        
        if (i < 18):
            npts = 2**18 - 2**i
            factor = 2**18 / 2**i
            yyi_pad = np.hstack((np.zeros((npts/2)), yyi, np.zeros((npts/2.))))
        else:
            factor = 1
            yyi_pad = yyi

        n_pad = len(yyi_pad)
        yyi_pad_shift = np.roll(yyi_pad, n_pad/2)
        new_spec = np.fft.ifft(yyi_pad_shift) #This is the smoothed spectra set for interpolation
        new_dv = delx/float(factor)
        
        new_x = np.arange(0,len(new_spec)/2, dtype=np.float64) * new_dv
        new_y = factor * new_spec[0: (len(new_spec)/2)] #removed -1 because of how python slices wrt IDL
        new_delx = np.mean(np.diff(new_x))
    else:
        new_x = x
        new_y = y
        new_delx = delx
    
    if AERI_dv / new_delx > 256:
        sfac = 256
    elif AERI_dv / new_delx > 128:
        sfac = 128
    elif AERI_dv / new_delx > 64:
        sfac = 64
    elif AERI_dv / new_delx > 32:
        sfac = 32
    elif AERI_dv / new_delx > 16:
        sfac = 16
    else:
        #If this happens, warn the user and still use sfac = 16
        sfac = 16
        print("WARNING in convolve2aeri: Unanticipated problem in computing new_x, setting sfac = 16")

    new_aeri_dv = AERI_dv / sfac
    max_v = np.max(new_x)
    new_aeri_wnum = np.arange(0,max_v / new_aeri_dv, 1, dtype=np.float64) * new_aeri_dv
    new_aeri_spec = np.interp(new_aeri_wnum, new_x, np.asarray(new_y, np.float64))

    #In our desire to have a spectrum with 2^n + 1 points, we may need to
    #throw away a few points (but these are probably in the taper)
    for i in np.arange(0,101):
        if 2**i >= len(new_aeri_wnum):
            break
    npts = 2**(i-1)
    new_aeri_wnum = new_aeri_wnum[0:npts+1] #add +1 b/c of how Python slices wrt IDL
    new_aeri_spec = new_aeri_spec[0:npts+1] #same as line before

    #Now fold this spectrum, and compute its interferogram
    n_aeri = len(new_aeri_spec)
    new_aeri_spec_folded = np.hstack((new_aeri_spec, new_aeri_spec[1:n_aeri - 1][::-1])) #Change to -1 b/c Py and IDL slicing
    n_fold = len(new_aeri_spec_folded)
    new_aeri_inter = np.fft.fft(new_aeri_spec_folded)
    new_aeri_inter = np.asarray(new_aeri_inter, dtype=np.float64)
    new_aeri_inter = np.roll(new_aeri_inter, -1*n_fold/2)
    new_aeri_opd = 1. / (2 * new_aeri_dv)
    new_aeri_xx = (np.arange(0, n_fold, 1, dtype=np.float64)/float(n_fold)*2 - 1) * new_aeri_opd
    
    #Now chop this at the desired optical path delay
    chop_idx = np.where((-1*AERI_OPD <= new_aeri_xx) & (new_aeri_xx < AERI_OPD))[0]
    nchop_idx = len(chop_idx)
    aeri_chop_inter = new_aeri_inter[chop_idx]

    #Then transform back into the spectral domain
    #print aeri_chop_inter
    n_chop = len(aeri_chop_inter)
    aeri_chop_inter_shift = np.roll(aeri_chop_inter, n_chop/2)
    final_aeri_spec = np.fft.ifft(aeri_chop_inter_shift)

    #Compute the scale factor that will account for the energy redistribution
    factor =len(final_aeri_spec) / float(len(new_aeri_inter))
    final_aeri_spec = factor * final_aeri_spec

    #And compute the wavenumber grid for this data
    final_aeri_dv = 1. / (2. * AERI_OPD)
    final_aeri_wnum = np.arange(0, len(final_aeri_spec)/2) * final_aeri_dv
    final_aeri_spec = final_aeri_spec[0:len(final_aeri_wnum)]

    #Final step: cut off data before and after the actual minimum and maximum wavenumber intervals of input data
    idxs = np.where((minv + tapersize <= final_aeri_wnum) & (final_aeri_wnum <= maxv-tapersize))[0]

    if len(idxs) <= 0 :
        print("ERROR")
        sys.exit()
    else:
        final_aeri_wnum = final_aeri_wnum[idxs]
        final_aeri_spec = final_aeri_spec[idxs]
    #dict = {'final_aeri_wnum': final_aeri_wnum, 'final_aeri_spec': final_aeri_spec}
    #np.savez('pyconv.npz', **dict)
    #apo_spec = apodizer.apodizer(final_aeri_spec)
    #n = len(final_aeri_spec)
   # imd = n/2
    #beer = np.zeros(n)
    #beer[0] = 1.
    #for i in np.arange(1,n/2+1):
    #    if i <= imd:
    #        beer[i] = (1 - ((i-1)/double(imd))**2)**2
    #    else:
    #        beer[i] = 0.
    
    #if (n % 2) == 0:
    #    beer[n/2:n] = beer[0:n/2][::-1]
    #else:
    #    beer[n/2:n] = beer[0: (n/2) + 1][::-1]

    #new_spec = np.fft.ifft(np.fft.fft(final_aeri_spec)*beer)

    return final_aeri_wnum, final_aeri_spec

def test():
    return
#if __name__ == '__main__':
    #ods = panel_file.panel_file('../idl_ex/out/ODdeflt_001', do_load_data=True)    
    #aw, aspec = convolve_to_aeri(ods.v, ods.data1)
    #idl = pidly.IDL()
    #conv = idl.convolve_to_aeri(ods.v, ods.data1)
    #np.savez('IDL_convolve.npz', **conv)
    #stop
    #plot(aw, aspec*10.)
    #show()
    #conv = np.load('ODdeflt_001.npz')
    #plot(ods.v, ods.data1)
    #plot(conv['wnum'], conv['spec'])
    #xlim([200,400])
    #show()
    #print conv['spec']/np.asarray(aspec, dtype=np.float64)
    #print len(conv['wnum']), len(conv['spec'])
    #print len(aw), len(aspec)
    #diff = conv['wnum'] - aw
    #print min(diff)
    #print max(diff)
    #diff = conv['spec'] - aspec
    #print min(diff)
    #print max(diff)
