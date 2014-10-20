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
from IPython.parallel import Client
import subprocess

def smoothAERI(a, dp, base_wnum, ods):
    idxs = np.where((base_wnum > a - dp) & (base_wnum < a + dp))[0]
    #ans = gauss(ods, idxs)
    ans = np.mean(np.exp(ods[:,idxs]*-1.), axis=1)
    return -1. * np.log(ans)

def read_and_interpODs(OD_dir):
   files = glob.glob(OD_dir + '/ODdeflt_*')

   #This code loads in the highest OD layer into memory, the wnum grid
   #for this layer becomes the standard wavenumber grid for all other layers
   files = np.sort(files)
   fd = panel_file.panel_file(files[0], do_load_data=True)
   base_wnum = fd.v
   ods = np.empty((len(files), len(base_wnum)))
   ods[len(files)-1] = fd.data1

   #This code loads in the rest of the OD layers and interpolated them to the
   #wavenumber grid that is of the highest layer
   for f in np.arange(len(files)-2,-1,-1):
       fd = panel_file.panel_file(files[f], do_load_data=True)
       ods[f] = np.interp(base_wnum, fd.v, fd.data1)

   return ods, base_wnum

def computeJacobian(spectra1, spectra2, deltaX):
    return (spectra1 - spectra2)/deltaX

class LBLPkg:
    def __init__(self, lblOUTdir):
        self.lbl_datadir = lblOUTdir
        ods, wnum = read_and_interpODs(lblOUTdir)
        z, t = t7r.readTape(lblOUTdir + '/TAPE7')
        self.ods = ods #Read in ODs here
        self.temp = t #Temperature profile
        self.z = z #height profile
        self.q = 3 #wv mixing ratio profile
        self.base_wnum = wnum #base wnum for the ODs 
        self.aeri_wnums = np.load('/home/greg.blumberg/python_pkgs/aeri_wnumgrid.npy')
    
    def getLBLdir(self):
        return self.lbl_datadir

    def trueJacobian(self, pert):
        wnum, fx = self.radianceTrue()
        print self.getLBLdir() + "../OUT_TPERT/TAPE7"
        try:
            reg_z, pert_temp = t7r.readTape(self.getLBLdir()+"../OUT_TPERT/TAPE7")
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
		    #print np.min(reg_ods-temporary_ods), np.max(reg_ods-temporary_ods)
		    #If the AERI observation suggests clouds are being viewed cloud optical depths may be added in here.

		    #Calculate the true way of doing the AERI Radiances
            print "Computing perturbed AERI radiance w/o apodization."
            true_aeri_wnum, true_fxp = self.radianceTrue(self.base_wnum, temp_temp, temporary_ods)
            print "Calculating Jacobian for height: ", str(reg_z[i])
            aeri_jacobian[i] = computeJacobian(fx, true_fxp, pert)
            true_fxprime[i] = true_fxp

        return aeri_jacobian, true_fxprime

    def applyDOD(self, dp, ods=None):
        if ods is None:
            ods = self.ods
        aeri_ods = np.empty((len(self.aeri_wnums), ods.shape[1]))
        for l in range(len(self.aeri_wnums)):
            aeri_ods[:,l] = smoothAERI(aeri_wnums[l], dp[l], base_wnum, ods)
        return aeri_ods

    def deltaODJacobian(self, dopfile, pert):
        dodp = np.load(dopfile)
        dp = dodp['dps']
        wnum, fx = radianceDeltaOD(dp)
        try:
            reg_z, pert_temp = t7r.readTape("../OUT_T/TAPE7")
            #reg_z, pert_mxr = t7r.readTape("../OUT_Q/TAPE7")
        except:
            raise Exception("Perturbation TAPE7s not in the correct directories ../OUT_T/ & ../OUT_Q")
        aeri_jacobian = np.empty((len(reg_z), len(wnum)))
        true_fxprime = np.empty((len(reg_z), len(wnum)))

        levels = np.arange(0, len(reg_z))
        for t in range(len(levels)):
            i = levels[t]
            print "Modifying profile at height: " + str(reg_z[i]) + ' km'
            temporary_ods = 1.*self.reg_ods
            temp_temp = 1.*self.temp
            temp_temp[i] = pert_temp[i]
            if i == 0:
                fd = panel_file.panel_file(self.lbl_datadir + '/ODdeflt_001', do_load_data=True)
                temporary_ods[0] = np.interp(self.base_wnum, fd.v, fd.data1)
                print "Bottommost layer, changing only the bottommost OD layer and temp[0]."
            elif i+1 == len(reg_temp):
                print self.lbl_datadir + '/ODdeflt_'  + str((3-len(str(i+1))) * '0' + str(i))
                fd = panel_file.panel_file(self.lbl_datadir + '/ODdeflt_'  + str((3-len(str(i+1))) * '0' + str(i)), do_load_data=True)
                temporary_ods[len(reg_temp)-1-1] = np.interp(self.base_wnum, fd.v, fd.data1)
                print "Top of the atmosphere: changing only the topmost OD layer and temp."
            else:
                #Replace the optical depths for the layer below level i (makes it layer i).
                fdbo = panel_file.panel_file(self.lbl_datadir + '/ODdeflt_' + str((3-len(str(i))) * '0' + str(i)), do_load_data=True)
                temporary_ods[i-1] = np.interp(self.base_wnum, fdbo.v, fdbo.data1)
                #Replace the optical depths for the layer above level i (makes it layer i+1)
                fdup = panel_file.panel_file(self.lbl_datadir + '/ODdeflt_' + str((3-len(str(i+1))) * '0' + str(i+1)), do_load_data=True)
                temporary_ods[i] = np.interp(self.base_wnum, fdup.v, fdup.data1)
                #print np.min(reg_ods-temporary_ods), np.max(reg_ods-temporary_ods)
            #If the AERI observation suggests clouds are being viewed cloud optical depths may be added in here.

            #Calculate the true way of doing the AERI Radiances
            print "Computing perturbed AERI radiance w/o apodization."
            true_aeri_wnum, true_fxp = self.radianceDeltaOD(dp, temp_temp, temporary_ods)
            print "Calculating Jacobian for height: ", str(reg_z[i])
            aeri_jacobian[i] = computeJacobian(fx, true_fxp, pert)
            true_fxprime[i] = true_fxp

        return aeri_jacobian, true_fxprime

    def radianceTrue(self, wnums=None, temp=None, ods=None):
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

    def monoRadiance(self, wnums=None, temp=None, ods=None):
        if wnums is None:
            wnums = self.base_wnum
        if temp is None:
            temp = self.temp
        if ods is None:
            ods = self.ods
        rad = rxf.rt(wnums, temp, ods)
        return wnums, rad

    def radianceDeltaOD(self, dp, temp=None, ods=None):
        if temp is None:
            temp = self.temp
        if ods is None:
            ods = self.ods
        DOD_ODs = applyDOD(dp, self.ods)
        wnum, rad = rxf.rt(self.aeri_wnums, self.temp, DOD_ODs)
        """Return Radiance"""
        return wnum, rad
