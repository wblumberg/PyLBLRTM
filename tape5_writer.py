from netCDF4 import Dataset
import sys
import numpy as np
from pylab import *

#Clausius Claperyon
def cc(temp, L):
    R_v = 461.5 # J/kg
    e_s0 = 6.11 #mb
    T_0 = 273.15 #K
    e = e_s0 * np.exp((L/R_v)*((1./T_0) - (1./(temp+273.15))))
    return e

def rh2w(rh, pres, temp):
    L_s = 2.83 * 10**6 #J/kg
    e_0 = .06112 #mb
    L = 2.5 * 10.**6 #J/kg
    e_s0 = 6.11 #mb

    e_s = cc(temp, L)
    e_i = cc(temp, L_s)

    e = (rh/100.) * e_s
    rh_ice = (e/e_i)*100.
    rh_liq = (e/e_s)*100.

    ice_idxs = np.where(temp <= -15)[0]
    rh[ice_idxs] = rh_ice[ice_idxs]

    mxr = ((0.622 * e)/(pres - e))*1000.

    return mxr, rh

"""
    sys.argv[1] should be the path to an ARM netCDF4 sonde file.
"""
#data = Dataset(sys.argv[1])
#out_file = sys.argv[2]

def makeFile(TAPE5_line, temperature, pressure, altitude, mxr, out_file):

    #print "Making TAPE5 for LBLRTM using: ", sys.argv[1]
    print "This TAPE5 will be saved in: ", out_file

    #temperature = data.variables['tdry'][:]
    #pressure = data.variables['pres'][:]
    #altitude = data.variables['alt'][:]/1000.
    #altitude = altitude - altitude[0]
    #rh = data.variables['rh'][:]

    #mxr, newrh = rh2w(rh, pressure, temperature)

    species_vec = np.ones(7,)

    full_atm = np.zeros((temperature.shape[0], 3+7))
    full_atm[:,3] = mxr

    fid = open(out_file, 'w')

    def write(format, var):
        try:
            fid.write(format % var)
        except:
            fid.write(format % 0)

    fid.write(TAPE5_line +'\n')
    fid.write('         1         2         3         4         5         6         7         8         9         0\n')
    fid.write('123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789 123456789\n')
    fid.write('$ None\n')
    """
    Section 1:  Using the LBLRTM to create the needed input for LBLDIS

    The LBLRTM is a line-by-line radiative transfer model that has been
    developed by Atmospheric and Environmental Research, Inc.  (AER).  The code,
    as well as the needed inputs, documentation, and examples, are available
    and freely available from their website: http://rtweb.aer.com.  The discussion
    below assumes that the reader is familiar with the LBLRTM and its terminology.

    The LBLRTM execution is controlled by the specification of the flags in the
    TAPE5 input file.  To create the output from the LBLRTM that is needed by
    the LBLDIS, these control flags need these specified values:

          IATM   = 1
          IEMIT  = 0
          IMRG   = 1
          IOD    = 0
          IPUNCH = 1

    Setting IATM to 1 indicates that the model should read in the profile
    information that is specified in the TAPE5 (i.e., to use the LBLATM module).

    Setting the IEMIT flag to 0 tells the LBLRTM to perform optical depth
    calculations only.

    Setting the IMRG flag to 1 will cause the LBLRTM to output the optical depths
    for each layer in their own file.

    Setting the IOD flag to 0 allows the LBLRTM to output the optical depths for
    each layer at the default spectral resolution that was used by the LBLRTM when
    it actually calculated these optical depths.  (This is the most efficient way
    to store the optical depth data, but requires the the LBLDIS code perform
    interpolation and/or integration to get the optical depth at each layer at the
    desired wavelength.)

    Setting the IPUNCH flag to 1 will result in the LBLRTM creating the TAPE7
    "punch" file, which contains the atmosphere layering and structure that was
    used in the optical depth calculations.  The TAPE7 file is read in by the
    LBLDIS code so that the LBLDIS routine uses the same values as the LBLRTM did.

    An example of a TAPE5 which has the options set correctly is given below:

    <SOF>
    $ Example TAPE5 file for LBLRTM to create input needed for LBLDIS
    HI=1 F4=1 CN=1 AE=0 EM=0 SC=0 FI=0 PL=0 TS=0 AM=1 MG=1 LA=0 MS=0 XS=0   00   00
     800.000  1000.000                                        0.0002    0.001
      6    2    5    1    1    7    1                                      360.000
       0.000     3.000     0.000
       0.000     0.500     1.000     2.000     3.000
    -1.
    \%\%\%

    """

    """ TAPE5 RECORD 1.2 VARIABLES """
    #HI
    IHIRAC = 1 # 1 - Voigt Profile
    #F4
    ILBLF4 = 1 # 1 - line-by-line bound is 25 cm-1 for all layer pressures
    #CN
    ICNTNM = 1 # 1 - Use the continumn
    #AE
    IAERSL = 0 # 0 - no aerosols used (this is also connected to the commented section below that has no aerosols)
    #EM
    IEMIT = 0 # if 0 - perform optical depth calculations only
    #SC
    ISCAN = 0 # 0 - Scanning function
    #IF
    IFILTR = 0 # 0 - no
    #PL
    IPLOT = 0 # 0 - no
    #TS
    ITEST = 0 # 0 - no
    #AM
    IATM = 1 # 1 - LBLRTM should read in the profile information specified in the TAPE5
    #MG
    IMRG = 1 # 1 - This will output the optical depths for each layer in their own files.
    #LA
    ILAS = 0 # Flag for laser options
    #MS
    IOD = 0 # 0 - output the ODs for each layer at the default spectral resolution used by the LBLRTM
    #XS
    IXSECT = 0 # 0 - use no cross-sections for the LBLRTM

    MPTS = 0 # number of optical depth values printed for the beginning and ending of each panel as a result of conv. for prev. layer
    NPTS = 0 # num of values printed for beginning and ending of each panel as a result of merge of current layer w/ prev layers

    """ TAPE5 RECORD 1.3 VARIABLES """
    V1 = 400 # Beginning wavenumber for the spectra
    V2 = 1900 # Ending wavenumber for the spectra
    SAMPLE = 4
    DVSET = 0 # THIS SELECTION MAY BE IMPORTANT FOR WHEN WE NEED TO GENERATE THE OPTICAL DEPTH DV CONSISTENT ACROSS ALL HEIGHTS
    ALFALO = 0.04 #DEFAULT = 0.04
    AVMASS = 36 # DEFAULT = 36
    DPTMIN = 0.0002
    DPTFAC = 0.001
    ILNFLG = 0 # DEFAULT = 0
    DVOUT = 0 #SELECTED DV GRID FOR THE OPTICAL DEPTH 'MONOCHROMATIC' OUTPUT SPACING (MUST BE LESS THAN OR EQ TO DEFAULT SPACING OF DVSET)
    NMOL_SCAL = 7 #NUMBER OF MOLECULES TO SCALE

    """ TAPE5 RECORD 1.3.A VARIABLES """
    #REQUIRED SINCE NMOL_SCAL = 7 IN THE ABOVE CODE
    HMOL_SCALE = '1M11111'

    #Molecules
    #1 - H2O ; water vapor
    #2 - CO2 ; carbon dioxide
    #3 - O3 ; ozone
    #4 - N2O ; nitrous oxide
    #5 - CO ; carbon monoxide
    #6 - CH4 ; methane
    #7 - O2 ; oxygen
    """ TAPE5 RECORD 1.3.B VARIABLES """
    #REQUIRED SINCE NMOL_SCAL = 7
    #Make sure you convert it to parts per million (not 380 but 380*10**-6
    HMOL_VALS = [1,380e-6,1,1,1,1,1]

    """ TAPE5 RECORD 3.1 VARIABLES """
    MODEL = 0 # Selects atmospheric profile
    ITYPE = 2 #selects type of path (slant path from H1 to H2)
    IBMAX = 49 # number of layer boundaries read in on Record 3.3B
    NOZERO = 1
    NOPRNT = 1
    NMOL = 7 # 7 molecules (max is 35)
    IPUNCH = 1 # 1 - will results in the LBLRTM creating TAPE7 which contains the layering and structure used in the OD calcs.
    #CO2MX = 380 #ppm (default is 330 ppm)

    """ TAPE5 RECORD 3.2 VARIABLES """
    H1 = 0 # Observer altitude
    H2 = 20 # End point altitude 
    ANGLE = 0 # zenith angle at H1 (degrees)

    """ TAPE5 RECORD 3.3B VARIABLES """
    ZNBD = [     0.000  ,  0.100   ,  0.200   ,  0.300   ,  0.400   ,  0.500  ,   0.600 ,    0.700,
         0.800  ,   0.900  ,   1.000   ,  1.250   ,  1.500  ,   1.750   ,  2.000 ,    2.250,
         2.500   ,  2.750  ,   3.000  ,   3.250   ,  3.500  ,   4.000   ,  4.500   ,  5.000,
         5.500 ,    6.000   ,  6.500  ,   7.000  ,   7.500  ,   8.000   ,  8.500   ,  9.000,
         9.500  ,  10.000  ,  10.500  ,  11.000  ,  11.500 ,   12.000  ,  12.500  ,  13.000,
        13.500  ,  14.000  ,  14.500   , 15.000   , 16.000  ,  17.000  ,  18.000  ,  19.000,
        20.000] #This is the array of heights in km the profile will be interpolated to in the LBLRTM
    
    ZNBD = ZNBD[::2]
    IBMAX = len(ZNBD)

    """ TAPE5 RECORD 3.5 & 3.6 VARIABLES """
    HMOD = " User supplied profile"
    JCHARP = "A" #Character representing Units for Pressure
    JCHART = "B" #Character representing units for temperature
    JLONG = ""
    JCHAR = " C666666"
    #altitude = []
    #temperature = []
    #pressure = []

    #species_vec = []

    #Card 1.2
    # Print out the main flags for the first line of the TAPE5 file
    write(' HI=%1i', IHIRAC)
    write(' F4=%1i', ILBLF4)
    write(' CN=%1i', ICNTNM)
    write(' AE=%1i', IAERSL)
    write(' EM=%1i', IEMIT)
    write(' SC=%1i', ISCAN)
    write(' FI=%1i', IFILTR)
    write(' PL=%1i', IPLOT)
    write(' TS=%1i', ITEST)
    write(' AM=%1i', IATM)
    write(' MG=%1i', IMRG)
    write(' LA=%1i', ILAS)
    write(' MS=%1i', IOD)
    write(' XS=%1i', IXSECT)
    write('  %2i', MPTS)
    write('  %2i', NPTS)
    fid.write('\n')

    #Card 1.3
    # This part is the wavenumber range (v1-v2) for the LBLRTM to operate with
    write('%10.3f', V1)
    write('%10.3f', V2)
    write('%10.3f', SAMPLE)
    write('%10.3f', DVSET)#DVSET
    write('%10.3f', ALFALO)
    write('%10.3f', AVMASS)
    write('%10.3f', DPTMIN)
    write('%10.3f', DPTFAC)
    write('    %1i', ILNFLG)
    write('     %10.3e', DVOUT)
    write('   %2i', NMOL_SCAL)

    fid.write('\n')

    #Card 1.3.a
    write('%s', HMOL_SCALE)
    fid.write('\n')

    #Card 1.3.b
    for i in range(1,len(HMOL_VALS)+1,1):
            write('%15.7e', HMOL_VALS[i-1])
            if i < len(HMOL_VALS) and round(i/8.) == i/8.:
                fid.write('\n')
    fid.write('\n')

    if IEMIT > 0:
    #Card 1.4
        write('%10.3f', TBOUND)
        write('%10.3f', SREMIS[1])
        fid.write('          ')
        fid.write('          ')
        write('%10.3f', SRREFL[1])
        fid.write('          ')
        fid.write('          ')
        fid.write('    ')
        write('%s', surf_refl)
        fid.write('\n')

    #Card 3.1
    write('%5i', MODEL)
    write('%5i', ITYPE)
    write('%5i', IBMAX)
    write('%5i', NOZERO)
    write('%5i', NOPRNT)
    write('%5i', NMOL)
    write('%5i', IPUNCH)
    #write('%2i', IFXTYP)
    #fid.write(' ')
    #write('%2i', MUNITS)
    #write('%10.3f', RE)
    #write('%10.3f', HSPACE)
    #write('%10.3f', VBAR)
    #write('%10.3f', CO2MX)
    #write('%10.3f', REF_L_AT)
    fid.write('\n')


    #Card 3.2
    # This is the range of heights (h1-h2) in the LBLRTM height grid 
    # and the angle the reciever is at 
    write('%10.3f', H1)
    write('%10.3f', H2)
    write('%10.3f', ANGLE)
    fid.write('\n')

    if abs(IBMAX) > 0:
        #Card 3.3B
        # 
        # This writes out the heights the user defined profile will be interpolated to by the LBLRTM
        for i in range(1,len(ZNBD)+1,1):
            write('%10.3f', ZNBD[i-1])
            if (i < len(ZNBD) and (round(i/8.) == i/8.)):
                fid.write('\n')
        fid.write('\n')
    else:
        #Record 3.3A
        """ we will not be using this section for the TAPE5s for the Fast Jacobian """
        AVTRAT = 0
        TDIFF1 = 0
        TDIFF2 = 0
        ALTD1 = 0
        ALTD2 = 0
        write('%10.3f', AVTRAT)
        write('%10.3f', TDIFF1)
        write('%10.3f', TDIFF2)
        write('%10.3f', ALTD1)
        write('%10.3f', ALTD2)
        fid.write('\n')
    #print MODEL


    if MODEL == 0:
        #Card 3.4
        # Printing out the number of record in the user-supplied profile
        IMMAX = len(altitude) #Gathering the number of records

        write('%5.0f', IMMAX)
        write('%s', HMOD) # HMOD is a 24 Character description of the profile
        fid.write('\n')
        #print IMMAX
        #Card 3.5 and Card 3.6
        for i in range(0,len(altitude),1):
            #Print out the height, pressure, and temperature values, along with the units string
            #print "Print"
            write('%10.3f', altitude[i])
            write('%10.3f', pressure[i])
            write('%10.3f', temperature[i])
            fid.write('     ')
            #Here is where we print out the units string
            write('%s', JCHARP)
            write('%s', JCHART)
            fid.write(' ')
            write('%s', JLONG)
            fid.write(' ')
            write('%s', JCHAR)
            fid.write('\n')

            #Print out the values for each atmospheric molecule constituent
            #print len(species_vec)
            #print len(species_vec)
            for j in range(1, len(species_vec)+1):
                if species_vec[j-1] == 1:

                    #Write out the value for atmospheric consituent j at height i
                    write('%10.3E', full_atm[i, j + 3 - 1])
                else:
                    fid.write('          ')
                if  (j < len(species_vec)) and (round(j/8.) == j/8.):
                    fid.write('\n')
            fid.write('\n')

    """
    if IAERSL > 0:
        #Card 4.1
        # Adding aerosol parameters for the IAERSL flag
        fid.write('%5i' % IHAZE)
        fid.write('%5i' % ISEASN)
        fid.write('%5i' % IVULCN)
        fid.write('%5i' % ICSTL)
        fid.write('%5i' % ICLD)
        fid.write('%5i' % IVSA)
        fid.write('          ')    #VIS
        fid.write('%10.3f' % WSS)
        fid.write('%10.3f' % WHH)
        fid.write('%10.3f' % RAINRT)
        fid.write('%10.3f' % GNDALT)
        fid.write('\n')

    if IXSECT > 0:    #!!! still needs to be developed
        # add in the x-section info
        fid.write('%5i%5i%5i selected x-sections are :\n' % (3, 0, 0))
        fid.write('CCL4      F11       F12 \\n')
        fid.write('%5i%5i  \\n' % [2, 0])
        fid.write('%10.3f     AAA\\n' % min(altitude))
        fid.write('%10.3e%10.3e%10.3e\\n' % [1.105e-04, 2.783e-09, 5.027e-04])
        fid.write('%10.3f     AAA\\n' % max(altitude))
        fid.write('%10.3e%10.3e%10.3e\\n' % [1.105e-04, 2.783e-09, 5.027e-04])


    if ISCAN == 1:
        #Card 8.1
        if ILS.ILS_use == 0:
            #Defaults
            HWHM = 0.01
            JV1 = V1
            JV2 = V2
            JEMIT = 1
            JFN = 0
            JVAR = 0
            SAMPL = 0
        else:
            HWHM = ILS.HWHM
            JV1 = ILS.V1
            JV2 = ILS.V2
            if get_flag(Flag_Vec, Flag_List, 'Rad_or_Trans') == 0:
                JEMIT = 0
            #JEMIT = 1;
            JFN = ILS.JFN
            JVAR = 0
            SAMPL = 0
        IUNIT = 12
        IFILST = 1
        NIFILS = 1
        JUNIT = 11

        fid.write('%10.6f' % HWHM)
        fid.write('%10.5f' % JV1)
        fid.write('%10.5f' % JV2)
        fid.write('   %2i' % JEMIT)
        fid.write('   %2i' % JFN)
        fid.write('   %2i' % JVAR)
        fid.write('%10.4f' % SAMPL)
        fid.write('   ')    #3X
        fid.write('%2i' % IUNIT)
        fid.write('   ')    #3X
        fid.write('%2i' % IFILST)
        fid.write('   ')    #3X
        fid.write('%2i' % NIFILS)
        fid.write('   ')    #3X
        fid.write('%2i' % JUNIT)
        fid.write('\n')
        fid.write('-1.\n')
        #end of Card 8.1

        blah = 0
        if blah == 1:
            #Reinterpolate scanned data to regular coordinate system
            #Card 1.1
            fid.write('$ Interpolation of scanned results\\n')

            #Card 1.2
            fid.write('%5i' % 0)        #IHIRAC
            fid.write('%5i' % 0)        #ILBLF4
            fid.write('%5i' % 0)        #ICNTNM
            fid.write('%5i' % 0)        #IAERSL
            fid.write('%5i' % 0)        #IEMIT
            fid.write('%5i' % 2)        #ISCAN
            fid.write('%5i' % 0)        #IFILTR
            fid.write('%5i' % 0)        #IPLOT
            fid.write('%5i' % 0)        #ITEST
            fid.write('%5i' % 0)        #IATM
            fid.write('%5i' % 0)        #IMRG
            fid.write('%5i' % 0)        #ILAS
            fid.write('%5i' % 0)        #IOD
            fid.write('%5i' % 0)        #IXSECT
            fid.write('%5i' % 0)        #MPTS
            fid.write('%5i' % 0)        #NPTS
            fid.write('\n')

            #Card 9.1
            fid.write('%10.6f' % HWHM)
            fid.write('%10.5f' % JV1)
            fid.write('%10.5f' % JV2)
            fid.write('   %2i' % 1)        #JEMIT
            fid.write('   %2i' % 1)        #JFN
            fid.write('   %2i' % 0)        #JVAR
            fid.write('          ')        #SAMPL
            fid.write('   ')        #3X
            fid.write('%2i' % 13)        #IUNIT
            fid.write('   ')        #3X
            fid.write('%2i' % IFILST)
            fid.write('   ')        #3X
            fid.write('%2i' % NIFILS)
            fid.write('   ')        #3X
            fid.write('%2i' % 11)        #JUNIT
            fid.write('\n')
            fid.write('-1.\n')
            #end of Card 9.1  
    """

    fid.write("%%%")
    print "Done."
