import numpy as np
import glob
import radxfer as rxf
import convolve2aeri as c2a
import sys
import panel_file
sys.path.append('../')
import apodizer
import tape7_reader as t7r
from pylab import *
#from IPython.parallel import Client
import subprocess
#from scipy import integrate

"""
    Object for Reading in LBLRTM data and doing performing calculations from the data.

    Last time this was run: (October 2013?)

    Written by: Greg Blumberg (OU/CIMMS)
    Email: wblumberg@ou.edu, greg.blumberg@noaa.gov

    TODO: Allow for upwelling and downwelling monochromatic calculations.
          Support calculating total flux calculations
          Perform convolution with boxcar filter function.
          Support the three methods of reducing the RTM calculations:
            1.) Convolve the raw optical depths before RTE
            2.) Transform ODs to Transmission then convolve then back to ODs then RTE
            3.) Compute layer to instrument transmission then convolve then back to ODs.
"""

def smoothAERI(a, dp, base_wnum, ods):
    """
        smoothAERI()

        smoothes the AERI spectra
    """
    idxs = np.where((base_wnum > a - dp) & (base_wnum < a + dp))[0]
    ans = np.mean(np.exp(ods[:,idxs]*-1.), axis=1)
    return -1. * np.log(ans)

def read_and_interpODs(OD_dir):
    """
        read_and_interpODs()

        Takes in the output directory from lblrun that has the ODdeflt_*
        files and reads in all of the optical depth files.  Because
        the LBLRTM TAPE5 writer writes the ODdeflt files out to different
        wavenumber grids, this function interpolates all of those to the 
        file with the maximum resolution.  This function returns
        the 2D array of optical depths and the maximum resolution wavenumber
        grid.
    """

    files = np.sort(glob.glob(OD_dir + '/ODdeflt_*'))[::-1]

    print "Reading in optical depth files..."
    #This code loads in the highest OD layer into memory, the wnum grid
    #for this layer becomes the standard wavenumber grid for all other layers
    fd = panel_file.panel_file(files[0], do_load_data=True)
    base_wnum = fd.v
    print "\tReading in:",files[0]
    ods = np.empty((len(files), len(base_wnum)))
    ods[len(files)-1] = fd.data1

    #This code loads in the rest of the OD layers and interpolated them to the
    #wavenumber grid that is of the highest layer
    for f in np.arange(1,len(files),1):
        print "\tReading in:",files[f]
        fd = panel_file.panel_file(files[f], do_load_data=True)
        ods[len(files)-1-f] = np.interp(base_wnum, fd.v, fd.data1)
    
    return ods, base_wnum

def computeJacobian(spectra1, spectra2, deltaX):
    return (spectra1 - spectra2)/deltaX

class LBLPkg:
    def __init__(self, lblOUTdir):
        self.lbl_datadir = lblOUTdir
        ods, wnum = read_and_interpODs(lblOUTdir)
        print "Reading in TAPE7..."
        z, t = t7r.readTape(lblOUTdir + '/TAPE7')
        print "LBLOUT files loaded."
        self.ods = ods #Read in ODs here
        self.temp = t #Temperature profile
        self.z = z #height profile
        self.q = 3 #wv mixing ratio profile interpolated from the LBLRTM TAPE7 (not implemented yet, don't know how to do this).
        self.base_wnum = wnum #base wnum for the ODs 
        #self.aeri_wnums = np.load('/home/greg.blumberg/python_pkgs/aeri_wnumgrid.npy')
    
    def getLBLdir(self):
        # Returns the LBLRUN output data directory
        return self.lbl_datadir

    def trueJacobian(self, pert):
        # Calculates the Jacobian  
        wnum, fx = self.radianceTrue()
        print self.getLBLdir() + "../OUT_TPERT/TAPE7"
        try:
            reg_z, pert_temp = t7r.readTape(self.getLBLdir()+"../OUT_TPERT/TAPE7")
            # The below is commented because my TAPE7 reader doesn't read in the interpolated WVMR grid yet.
            #reg_z, pert_mxr = t7r.readTape(self.getLBLdir()+"../OUT_QPERT/TAPE7")
        except:
            raise Exception("Perturbation TAPE7s not in the correct directories ../OUT_T/ & ../OUT_Q")
        aeri_jacobian = np.empty((len(reg_z), len(wnum)))
        true_fxprime = np.empty((len(reg_z), len(wnum)))

        levels = np.arange(0, len(reg_z))
        for t in range(len(levels)):
            i = levels[t]
            print "Modifying profile at height: " + str(reg_z[i]) + ' km'
            temporary_ods = 1.*self.ods
            temp_temp = 1.*self.temp
            temp_temp[i] = pert_temp[i]

            if i == 0:
                fd = panel_file.panel_file(self.lbl_datadir + '/ODdeflt_001', do_load_data=True)
                temporary_ods[0] = np.interp(self.base_wnum, fd.v, fd.data1)
                print "Bottommost layer, changing only the bottommost OD layer and temp[0]."
            elif i+1 == len(self.temp):
                print self.lbl_datadir + '/ODdeflt_'  + str((3-len(str(i+1))) * '0' + str(i))
                fd = panel_file.panel_file(self.lbl_datadir + '/ODdeflt_'  + str((3-len(str(i+1))) * '0' + str(i)), do_load_data=True)
                temporary_ods[len(self.temp)-1-1] = np.interp(self.base_wnum, fd.v, fd.data1)
                print "Top of the atmosphere: changing only the topmost OD layer and temp."
            else:
            #Replace the optical depths for the layer below level i (makes it layer i).
                fdbo = panel_file.panel_file(self.lbl_datadir + '/ODdeflt_' + str((3-len(str(i))) * '0' + str(i)), do_load_data=True)
                temporary_ods[i-1] = np.interp(self.base_wnum, fdbo.v, fdbo.data1)
            #Replace the optical depths for the layer above level i (makes it layer i+1)
                fdup = panel_file.panel_file(self.lbl_datadir + '/ODdeflt_' + str((3-len(str(i+1))) * '0' + str(i+1)), do_load_data=True)
                temporary_ods[i] = np.interp(self.base_wnum, fdup.v, fdup.data1)
		    #If the AERI observation suggests clouds are being viewed cloud optical depths may be added in here.

		    #Calculate the true way of doing the AERI Radiances
            print "Computing perturbed AERI radiance w/o apodization."
            true_aeri_wnum, true_fxp = self.radianceTrue(self.base_wnum, temp_temp, temporary_ods)
            print "Calculating Jacobian for height: ", str(reg_z[i])
            aeri_jacobian[i] = computeJacobian(fx, true_fxp, pert)
            true_fxprime[i] = true_fxp

        return aeri_jacobian, true_fxprime
    
    def aeri_radiance(self, wnums=None, temp=None, ods=None):
        """
            aeri_radiance()

            Calculates the radiance values that would be observed by the AERI using the LBLRTM output data.
            Unless the arguments above are specified, this function uses the optical depth and temperature data
            loaded in through this LBLPkg object.
        """
        if wnums is None:
            wnums = self.base_wnum
        if temp is None:
            temp = self.temp
        if ods is None:
            ods = self.ods
        rad = rxf.rt(wnums, temp, ods)
        wnum, rad = c2a.convolve_to_aeri(wnums, rad)
        rad = apodizer.apodizer(rad)
        idxs = np.where((wnum >= self.aeri_wnums[0]) & (wnum <= (self.aeri_wnums[-1] + .1)))[0]
        rad = rad[idxs]
        return wnum[idxs], rad

    def monoRadiance(self, zenith_angle=0, sfc_t=None, sfc_e=None, upwelling=False):
        """
            monoRadiance()

            Calculates the monochromatic radiance depending upon certain parameters.
            By default, it calculates the downwelling monochromatic radiance values.

            Parameters
            ----------
            zenith_angle : zenith angle of the calculation (degrees; default=0)
            sfc_t : surface temperature (Kelvin; default=None)
            sfc_e : surface emissivity (unitless)
            upwelling : switch to compute upwelling vs downwelling radiance
                        False - downwelling
                        True - upwelling
        """
        wnums = self.base_wnum
        temp = self.temp
        ods = self.ods
        rad = rxf.rt(wnums, temp, ods, zenith_angle=zenith_angle, sfc_t=sfc_t, sfc_e=sfc_e, upwelling=upwelling)
        return wnums, rad

    def od2trans(self, od=None):
        if od is None:
            od = self.ods
        return np.exp(-od)

    def computeFlux(self, upwelling=False, v1=400, v2=2000, sfc_t=300, sfc_e=1, dtheta=30):
        # if layer = 0
        #   compute TOA flux
        # if layer = 1
        #   compute BOA flux
        #
        # integrate I_v(theta) sin(theta) cos(theta) d(theta) d(phi) d(wavenumber)
            #wnum = self.base_wnum
        print "Upwelling?:", upwelling
        range_of_thetas = np.arange(0,90 + dtheta,dtheta)
        radiances = np.empty((len(range_of_thetas), len(self.base_wnum)))
        for i in range(len(range_of_thetas)):
            theta = range_of_thetas[i]
            #print theta
            #print self.monoRadiance(theta, sfc_t, sfc_e, upwelling)
            #radiances[i,:] = self.monoRadiance(theta, sfc_t, sfc_e, upwelling)[0]# * np.sin(np.radians(theta)) * np.cos(np.radians(theta))
            #plot(self.base_wnum, radiances[i,:])
            radiances[i,:] = self.monoRadiance(theta, sfc_t, sfc_e, upwelling)[1] * np.sin(np.radians(theta)) * np.cos(np.radians(theta))
            #plot(self.base_wnum, radiances[i,:])
        show()
        print radiances.shape
        # After integrating over theta
        integrated = np.trapz(radiances, range_of_thetas, dx=dtheta, axis=0)
        print integrated
        # After integrating over phi
        integrated = integrated * (2*np.pi) 
        print integrated 
        print integrated.shape 
        #print integrate.quad(lambda x: np.interp(x, self.base_wnum, integrated), v1, v2)[0] * 0.001

