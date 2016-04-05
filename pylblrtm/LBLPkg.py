import numpy as np
import glob
import radxfer as rxf
import convolve2aeri as c2a
import sys
import panel_file
sys.path.append('../')
import apodizer
import tape7_reader as t7r
import subprocess
from scipy import convolve

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
    switch = 1

    files = np.sort(glob.glob(OD_dir + '/ODdeflt_*'))[::-1]

    print "Reading in optical depth files..."
    #This code loads in the highest OD layer into memory, the wnum grid
    #for this layer becomes the standard wavenumber grid for all other layers
    if switch == 0:
        fd = panel_file.panel_file(files[0], do_load_data=True)
        base_wnum = fd.v
        print "\tReading in:",files[0]
        ods = np.empty((len(files), len(base_wnum)))
        ods[len(files)-1] = fd.data1
        begin_idx = 1
    else:
        fd = panel_file.panel_file(files[-1], do_load_data=True)
        base_wnum = fd.v
        ods = np.empty((len(files), len(base_wnum)))
        begin_idx = 0

    #This code loads in the rest of the OD layers and interpolated them to the
    #wavenumber grid that is of the highest layer
    for f in np.arange(begin_idx,len(files),1):
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

    def filterSpectra(self, y, wnum, delta=5):
        """
            filterSpectra

            Convolves a spectrum with the boxcar function specified by the delta
            argument.  Returns the spectrum back to the user once done.

            Arguments
            ---------
            y - the spectrum requested (on a specific wavenumber grid)
            wnum - the evenly spaced wavenumber grid (to determine the window of the boxcar)
            delta - the size of the boxcar window (cm-1)
        """
        dv = np.diff(wnum)[0]
        N = delta/dv
        r = convolve(y, np.ones(N)/N, mode='valid')
        return r

    def filterRadiance(self, delta=5, method=1, zenith_angle=0, sfc_t=None, sfc_e=None, upwelling=False, debug=False):
            
        wnums = self.base_wnum
        filtered_wnum = self.filterSpectra(wnums, wnums, delta)
        filtered_ods = np.empty((len(self.ods), len(filtered_wnum)))
        if method == 1:
            # Do convolution in OD space.
            for i in range(len(self.ods)):
                if debug is True:
                    print "Convolving ODs @ LEVEL:", i
                    print "Max:", np.max(self.ods[i,:]), "Min:", np.min(self.ods[i,:])
                filtered_ods[i,:] = self.filterSpectra(self.ods[i,:], wnums, delta)
        elif method == 2:
            # Do convolution in transmission space.
            for i in range(len(self.ods)):
                trans = self.od2t(self.ods[i,:])
                if debug is True:
                    print "Convolving transmission @ level:",i
                    print "Max:", np.max(trans), "Min:", np.min(trans)
                trans_cov = self.filterSpectra(trans, wnums, delta)
                filtered_ods[i,:] = self.t2od(trans_cov)
        elif method == 3:
            # Do the convolution in layer-to-instrument transmission.
            
            # Convert ODs to transmission
            trans = self.od2t(self.ods)
            if upwelling is True:
                trans = trans[::-1,:] # Reverse the transmissivity 
            new_l2i_t = np.empty((len(self.ods), len(filtered_wnum)))
            new_trans = np.empty((len(self.ods), len(filtered_wnum)))
            for i in xrange(len(new_trans)):
                # Compute layer-to-instrument transmission.
                l2i_t = np.prod(trans[:i,:], axis=0)
                # Convolve that mother.
                new_l2i_t[i,:] = self.filterSpectra(l2i_t, wnums, delta)  
            # Set the first level equal to the layer to instruemtn transmissivity
            new_trans[0,:] = new_l2i_t[0,:]
            for i in xrange(1, len(new_trans)):
                # Divide current layer L2I Transmissivity by the layer "below" that.
                # e.g. (trans_1 * trans_2 * trans_3) / (trans_1 * trans_2) = trans_3
                new_trans[i,:] = np.divide(new_l2i_t[i,:], new_l2i_t[i-1,:])

            # Convert back to the optical depth space.
            filtered_ods = self.t2od(new_trans)
            if upwelling is True:
                # If we're computing upwelling, reverse the optical depths again
                filtered_ods = filtered_ods[::-1,:] 
        else:
            print "Method not supported."
            return None, None      
        
        rad = rxf.rt(filtered_wnum, self.temp, filtered_ods, zenith_angle=zenith_angle, sfc_t=sfc_t, sfc_e=sfc_e, upwelling=upwelling)

        return filtered_wnum, rad

    def od2t(self, od=None):
        if od is None:
            od = self.ods
        return np.exp(-od)
    
    def t2od(self, t):
        return -np.log(t)

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

