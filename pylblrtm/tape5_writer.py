from netCDF4 import Dataset
import sys
import numpy as np

"""
    SCRIPT NAME:
    tape5_writer.py

    AUTHOR:
    Greg Blumberg (OU/CIMMS)
    wblumberg@ou.edu, greg.blumberg@ou.edu

    DESCRIPTION:
    This file will generate a default TAPE5 file using these arguments:
    
    TAPE5_line - a line describing what the TAPE5 file does
    temperature - an array of length l (Celsius), which can come from a radiosonde
    pressure - an array of length l (millibars)
    altitude - an array of length l (kilometers)
    mxr - an array of length l of the water vapor mixing ratio (g/kg) 
    outfile - the name of the file to write all the TAPE5 calculations to.

    Based off of Dave Turner's rundecker, but significantly more stripped down
"""

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

def makeFile(out_file, V1, V2, MODEL, ZNBD=None, IEMIT=0, HMOL_VALS=[1,380e-6,1,1,1,1,1], upwelling=False,**kwargs):
    """
        makeFile()

        This function creates an LBLRTM TAPE5 input file using the arguments passed to it.
        This TAPE5 writer only recognizes the contributions from the first 7 molecules of the 
        TAPE1 file, which are H2O, CO2, O3, N2O, CO, CH4, and O2.

        This TAPE5 writer is a work in progress.  Not all LBLRTM functions are supported, but
        currently this function can create TAPE5 files for LBLRTM outputs that can be read by
        LBLDIS or the LBLPkg object included in PyLBLRTM.  This function does this through
        the IEMIT=0 switch, which outputs gaseous optical depth files from the LBLRTM instead of
        radiance files.  In addition, switching IEMIT to 1 will switch the program to perform 
        it's own RT calculations and output them into TAPE10 and TAPE12 files.  It might be nice
        in the future to allow the LBLPkg to read in those files too.

        Although the LBLRTM supports scanning, interpolating, and FFTing the radiance spectrum
        it computes, this function does not support that capability yet.  At some point,
        I'll implement it and make a class that is an ILS (instrument line shape) class that
        contains all the parameters the LBLRTM needs to convert the monochromatic radiance to
        something that is observed by an instrument.

        Website explaining the TAPE5 format:
        http://web.gps.caltech.edu/~drf/misc/lblrtm/lblrtm_toc.html

        Required Arguments
        ------------------
        out_file : a string containing the filename and/or path for the TAPE5 file.
        V1 : the beginning wavenumber of the portion of the spectra we want to simulate (cm-1)
        V2 : the ending wavenumber of the portion of the spectra we want to simulate (cm-1)
        MODEL : an integer describing whether or not the TAPE5 should use a user-defined profile
                or a default profile that is managed by the LBLRTM
                
                Options:
                    0 - user-defined profile (requries that wvmr, pres, tmpc, and hght are specified)
                    1 - tropical model
                    2 - midlatitude summer model
                    3 - midlatitude winter model
                    4 - subarctic summer model
                    5 - subarctic winter model
                    6 - U.S. standard 1976 

        Optional Arguments With Defaults
        --------------------------------
        HMOL_VALS : the scaling factor for the 7 primary gases.  If the value is not 1,
                    the value is specified in parts per million * e-6.  Default values if not
                    specified is 380e-6 for CO2 and unscaled for all other gases.

        ZNBD : an array containing the height grid for the LBLRTM to use (km) (OPTIONAL)
               not setting this variable will cause the writer to use a default array going up to
               20 km.
        
        IEMIT : configure TAPE5 to perform optical depth calculations or radiance/transmittance calculations
                0 - optical depth only (default, used for LBLDIS or your own RT calculations)
                1 - radiance/transmittance (can apply instrument lineshape)
        
        upwelling : if IEMIT is equal to 1, computes the upwelling radiance if True
                Default is False
                True - computes upwelling radiance (sets H1 - max(ZNBD), H2 - min(ZBND), ANGLE=180)
                False - computes downwelling radiance (sets H1 - min(ZBND), H2 - max(ZNBD), ANGLE=0)
                if True, checks for sfc_e and sfc_t in kwargs.
                if False, sets sfc_e = 1 and sfc_1 = 0.0001
                 
        Keyword Arguments
        -----------------
        TAPE5_line : a string describing what the TAPE5 does or is for (for record keeping)
        IXSECT : use heavy molecule cross sections (0 - no, 1 - yes)
        ISCAN : convolve result with AERI ILS and output to TAPE27 (rad) and TAPE28 (trans)
                0: no, 3: yes 
             
        if MODEL == 0:
            wvmr : an array containing the water vapor mixing ratio profile (g/kg)
            pres : an array containing the pressure profile (mb)
            tmpc : an array containing the temperature profile (C)
            hght : an array containing the height profile (km)

        if IEMIT != 0:
            # must be specified if an upwelling calculation is done
            sfc_r : what type of surface reflectance to consider
                    's' - spectral
                    'l' - Lambertian (only applicable for upwelling calculations 90 <= angle <= 180)
            sfc_e : emissivity value at H2 (point furthest away from observer)
            sfc_t : surface temperature at H2 (Kelvin)
            ILS : an instrument lineshape object to apply to the spectrum (Not implemented yet)
    """

    print "This TAPE5 will be saved in: ", out_file

    species_vec = np.ones(7,)

    # Initialize the full_atm array to
    if MODEL == 0:
        if 'wvmr' not in kwargs or 'tmpc' not in kwargs or 'pres' not in kwargs or 'hght' not in kwargs:
            print "MODEL == 0, but argument is missing kwargs to specify the user-specified profile."
            print "TAPE5 not created."
            return

        temperature = kwargs.get('tmpc')
        full_atm = np.zeros((temperature.shape[0], 3+7))
        full_atm[:,3] = kwargs.get('wvmr') 
        pressure = kwargs.get('pres')
        altitude = kwargs.get('hght')

    # Open the up the file to write to.
    fid = open(out_file, 'w')

    # Create a function that will write to the TAPE5 file using a specific format 
    def write(format, var):
        try:
            fid.write(format % var)
        except:
            fid.write(format % 0)

    # Write the line describing the TAPE5 and what it is for
    fid.write(kwargs.get('TAPE5_line', '') +'\n')

    # Write some text file alignment text
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
    # This contains the primary settings for the LBLRTM.  It tells the model how we want to run the code.
    #HI
    IHIRAC = 1 # 1 - Voigt Profile
    #F4
    ILBLF4 = 1 # 1 - line-by-line bound is 25 cm-1 for all layer pressures
    #CN
    ICNTNM = 1 # 1 - Use the continumn
    #AE
    IAERSL = 0 # 0 - no aerosols used (this is also connected to the commented section below that has no aerosols)
               # aerosol calculations handled by LBLDIS
    #EM
    IEMIT = IEMIT # 0 - perform optical depth calculations only
                  # 1 - radiance and transmittance (W / cm2 sr-1 cm-1) 
    #SC
    ISCAN = kwargs.get('ISCAN', 0) # 0 - No scanning function
                                   # 1 - apply scanning function
                                   # 2 - apply interpolation procedure to results
                                   # 3 - apply FFT scan
    #IF
    IFILTR = 0 # 0 - no, 1 - yes
    #PL
    IPLOT = 0 # 0 - no, 1 - yes
    #TS
    ITEST = 0 # 0 - no, 1 - yes
    #AM
    IATM = 1 # 1 - LBLRTM should read in the profile information specified in the TAPE5
    #MG
    if IEMIT == 0: # When we want to look at the ODs and not the radiances/transmission.
        IMRG = 1 # 1 - This will output the optical depths for each layer in their own files.
    else:
        IMRG = 0 # This will ensure that the TAPE10 and TAPE12 will have something in it.
    #LA
    ILAS = 0 # Flag for laser options
    #MS
    IOD = 0 # 0 - output the ODs for each layer at the default spectral resolution used by the LBLRTM
            # 
    #XS
    IXSECT = kwargs.get('IXSECT', 0) # 0 - use no cross-sections for the LBLRTM, 1 - use cross-sections

    MPTS = 0 # number of optical depth values printed for the beginning and ending of each panel as a result of conv. for prev. layer
    NPTS = 0 # num of values printed for beginning and ending of each panel as a result of merge of current layer w/ prev layers

    """ TAPE5 RECORD 1.3 VARIABLES """
    # This describes the portion of the EM spectrum we want to look at.
    V1 = V1 # Beginning wavenumber for the spectra
    V2 = V2 # Ending wavenumber for the spectra
    if V1> V2:
        print 'INVALID WAVELENGTH RANGE.'
        return
    SAMPLE = 4
    DVSET = 0 # THIS SELECTION MAY BE IMPORTANT FOR WHEN WE NEED TO GENERATE THE OPTICAL DEPTH DV CONSISTENT ACROSS ALL HEIGHTS
    ALFALO = 0.04 #DEFAULT = 0.04
    AVMASS = 36 # DEFAULT = 36
    DPTMIN = -1
    DPTFAC = 0.001
    ILNFLG = 0 # DEFAULT = 0
    DVOUT = 0 #SELECTED DV GRID FOR THE OPTICAL DEPTH 'MONOCHROMATIC' OUTPUT SPACING (MUST BE LESS THAN OR EQ TO DEFAULT SPACING OF DVSET)
    NMOL_SCAL = 7 #NUMBER OF MOLECULES TO SCALE

    """ TAPE5 RECORD 1.3.A VARIABLES """
    #REQUIRED SINCE NMOL_SCAL = 7 IN THE ABOVE CODE
    #HMOL_SCALE = '1M11111' # Default
    HMOL_SCALE = ''
    for i in HMOL_VALS:
        if i == 1:
            # scaling factor used directly to scale profile
            HMOL_SCALE += '1'
        else:
            # volume mixing ratio wrt dry air for the total column to which the profile will be scaled
            HMOL_SCALE += 'M'

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
    #Make sure you convert it to parts per million (not 380 but 380*10**-6)..will scale CO2?
    #This is now set elsewhere but default will be: HMOL_VALS = [1,380e-6,1,1,1,1,1]

    """ TAPE5 RECORD 1.4 VARIABLES """
    # emissivity and t boundary and reflectivity variables (describes sfc of planet)
    if upwelling is False:
        # Set parameters for a downwelling calculation
        TBOUND = 0.00001 # Kelvin
        SREMIS = 1 # surface emissivity 
        surf_refl = 's' # treat the surfaces as Lambertian (reflects isotropically)
    else:
        # Set parameters for an upwelling calculation
        TBOUND = kwargs.get('sfc_t', None)
        SREMIS = kwargs.get('sfc_e', 1)
        if TBOUND is None:
            print "You asked for an upwelling calculation, but didn't give a SFC_T argument to the function!"
            print "Aborting TAPE5 creation."
            sys.exit()
        surf_refl = kwargs.get('sfc_r', 'l')
        
    """ TAPE5 RECORD 3.1 VARIABLES """
    # This describes the type of atmosphere, the path the radiation goes to, and the gases in the simulation.
    #MODEL = 0 # Selects atmospheric profile
    ITYPE = 2 #selects type of path (slant path from H1 to H2)
    IBMAX = 49 # number of layer boundaries read in on Record 3.3B
    NOZERO = 1
    NOPRNT = 1
    NMOL = 7 # 7 molecules (max is 35)
    IPUNCH = 1 # 1 - will results in the LBLRTM creating TAPE7 which contains the layering and structure used in the OD calcs.
    #CO2MX = 380 #ppm (default is 330 ppm)

    """ TAPE5 RECORD 3.3B VARIABLES """
    # This tells the LBLRTM what kind of height grid it should interpolate the T/P/gas concentration data to.
    if ZNBD is None:
        ZNBD = [     0.000  ,  0.100   ,  0.200   ,  0.300   ,  0.400   ,  0.500  ,   0.600 ,    0.700,
                     0.800  ,   0.900  ,   1.000   ,  1.250   ,  1.500  ,   1.750   ,  2.000 ,    2.250,
                     2.500   ,  2.750  ,   3.000  ,   3.250   ,  3.500  ,   4.000   ,  4.500   ,  5.000,
                     5.500 ,    6.000   ,  6.500  ,   7.000  ,   7.500  ,   8.000   ,  8.500   ,  9.000,
                     9.500  ,  10.000  ,  10.500  ,  11.000  ,  11.500 ,   12.000  ,  12.500  ,  13.000,
                    13.500  ,  14.000  ,  14.500   , 15.000   , 16.000  ,  17.000  ,  18.000  ,  19.000,
                    20.000] #This is the array of heights in km the profile will be interpolated to in the LBLRTM
        ZNBD = np.concatenate((np.arange(11)*0.1,np.arange(10)*0.25+1.25,\
                          np.arange(23)*0.5+4.0,np.arange(5)+16, np.arange(10)*2+22, np.arange(8)*4+42))
    """ TAPE5 RECORD 3.2 VARIABLES """
    # This tells the LBLRTM things about where the "observer" for the radiation is and what angle it is observing at.
    if upwelling is False:
        # set parameters for a downwelling calculation
        H1 = np.min(ZNBD) # Observer altitude
        H2 = np.max(ZNBD) # End point altitude 
        ANGLE = 0 # zenith angle at H1 (degrees) (can be used for satellite computations)
    else:
        # set parameters for an upwelling calculation
        H1 = np.max(ZNBD)
        H2 = 0.00001
        ANGLE = 180

    IBMAX = len(ZNBD)

    """ TAPE5 RECORD 3.5 & 3.6 VARIABLES """
    # This tells the LBLRTM things about the user-supplied profile (if one is given; MODEL == 0)
    HMOD = " User supplied profile"
    JCHARP = "A" #Character representing Units for Pressure (mb)
    JCHART = "B" #Character representing units for temperature (Celsius
    JLONG = ""
    print "For all other gases (except for WVMR), will be using the US Standard ATM profile."
    JCHAR = " C666666" # This represents the units of the gas concentrations we are specifying (C is g/kg), 6 means use the US Std Atmos.

    print "Now writing lines to the TAPE5."

    #Card 1.2
    print "Writing Card 1.2..."
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
    write(' XS=%1i', IXSECT) # This was some kind of cross sections
    write('  %2i', MPTS)
    write('  %2i', NPTS)
    fid.write('\n')

    #Card 1.3
    print "Writing Card 1.3..."
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
    print "Writing Card 1.3.a..."
    write('%s', HMOL_SCALE)
    fid.write('\n')

    #Card 1.3.b - Write out the molecular scaling values
    print "Writing Card 1.3.b..."
    for i in range(1,len(HMOL_VALS)+1,1):
            write('%15.7e', HMOL_VALS[i-1])
            if i < len(HMOL_VALS) and round(i/8.) == i/8.:
                fid.write('\n')
    fid.write('\n')

    if IEMIT > 0:
        print "Writing Card 1.4..."
        #Card 1.4...this card must be written if the LBLRTM is to do radiance/transmissivity calculations
        write('%10.3f', TBOUND)
        write('%10.3f', SREMIS)
        fid.write('          ')
        fid.write('          ')
        fid.write('          ')
        #write('%10.3f', SRREFL[1])
        fid.write('          ')
        fid.write('          ')
        fid.write('    ')
        write('%s', surf_refl)
        fid.write('\n')

    #Card 3.1
    print "Writing Card 3.1..."
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
    print "Writing Card 3.2..."
    # This writes the range of heights (h1-h2) in the LBLRTM height grid 
    # and the angle the reciever is at 
    write('%10.3f', H1)
    write('%10.3f', H2)
    write('%10.3f', ANGLE)
    write('%10.3f', np.max([H1, H2]))  # RANGE parameter
    fid.write('\n')

    if abs(IBMAX) > 0:
        print "Writing Card 3.3.b..."
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

    # This if statement only is entered if we've selected the user-defined profile.
    # If we're using other default atmospheres (e.g. US Standard Atmosphere) then this isn't run.
    # However, if the user-defined profile is defined, then all gasses except for H2O are assumed
    # to come from the US standard atmosphere profile (setting 6).
    if MODEL == 0:
        #Card 3.4
        print "User supplied a profile...writing Card 3.4..."
        # Printing out the number of record in the user-supplied profile
        
        idx = np.where(ZNBD > np.max(altitude))[0]
        IMMAX = len(idx) + len(altitude) #Gathering the number of records

        write('%5.0f', IMMAX)
        write('%s', HMOD) # HMOD is a 24 Character description of the profile
        fid.write('\n')
        #print IMMAX
        #Card 3.5 and Card 3.6
        for i in range(0,len(altitude),1):
            #Print out the height, pressure, and temperature values, along with the units string
            #print "Print"
            write('%10.3f', altitude[i]) # units are km
            write('%10.3f', pressure[i]) # units are mb
            write('%10.3f', temperature[i]) # units are Celsius
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
            for j in range(1, len(species_vec)+1):
                if species_vec[j-1] == 1:

                    #Write out the value for atmospheric consituent j at height i
                    write('%10.3E', full_atm[i, j + 3 - 1])
                else:
                    fid.write('          ')
                if  (j < len(species_vec)) and (round(j/8.) == j/8.):
                    fid.write('\n')
            fid.write('\n')
        # If the supplied profile doesn't go high enough, fill in the rest of the profile
        # using the US standard atmosphere.
        print "Supplied profile doesn't go high enough...filling in the profile with US Std Atmos..."
        if np.max(altitude) < np.max(ZNBD):
            idx = np.where(ZNBD > np.max(altitude))[0]
            for i in idx:
                #Print out the height, pressure, and temperature values, along with the units string
                #print "Print"
                write('%10.3f', ZNBD[i]) # units are km
                write('%10.3f', 0) # units are mb
                write('%10.3f', 0) # units are Celsius
                fid.write('     ')
                #Here is where we print out the units string
                write('%s', 6)
                write('%s', 6)
                fid.write(' ')
                write('%s', JLONG)
                fid.write(' ')
                write('%s', ' 6666666')
                fid.write('\n')

                #Print out the values for each atmospheric molecule constituent
                for j in range(1, len(species_vec)+1):
                    if species_vec[j-1] == 1:
                        #Write out the value for atmospheric consituent j at height i
                        write('%10.3E', 0)
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
    """
    if IXSECT > 0:    #!!! still needs to be developed
        # add in the x-section info
        # Got these x-section info from DDT's rundecker.pro on 8/26/2016
        fid.write('%5i%5i%5i selected x-sections are :\n' % (3, 0, 0))
        fid.write('CCL4      F11       F12 \n')
        fid.write('%5i%5i XS 1995 UNEP values\n' % (2, 0))
        fid.write('%10.3f     AAA\n' % min(altitude))
        fid.write('%10.3E%10.3E%10.3E\n' % (1.105e-04, 2.783e-04, 5.027e-04))
        fid.write('%10.3f     AAA\n' % max(altitude))
        fid.write('%10.3E%10.3E%10.3E\n' % (1.105e-04, 2.783e-04, 5.027e-04))
    ISCAN=3

    opd = 1.03702766 # AERI
    delv = 1./(2*opd)
    npts = long(10000/delv+1)
    varray = np.arange(npts, dtype=int)*delv + delv
    idx1 = np.where(varray > V1 + 50)[0]
    idx2 = np.where(varray > V2 - 50)[0]

    #ISCAN = 1; SCANFN, ISCAN = 2; INTRPL, ISCAN = 3; FFTSCN
    # Currently only set to do FFTSCN for AERI
    if ISCAN == 3:  # 3 means use FFT
        print "Writing Card 10.1..."
        print "WARNING!  This portion of the code hasn't been tested.  Use with caution."
        #Card 8.1
        #Defaults to using AERI ILS characteristics
        #HWHM = 1.03702766
        #JV1 = 450.32550
        #JV2 = 1750.19440
        #JEMIT = 1
        #JFNin = 1 # 1-triangle scanning function, 0-boxcar scanning function
        #MRATin = -4 # no prescanning boxcaring is performed 
        #DVOUT = 0.482
        #IUNIT = 12 # TAPE12 to be scanned
        #IFILST = 1
        #NIFILS = 1
        #JUNIT = 13 # Output file for the scanned results
        # Will write out to TAPE13 (radiance) and TAPE14 (transmittance)
        # with the AERI ILS convolved. Settings taken from DDT rundecker.pro
        string = '%10.8f%10.5f' % (opd, varray[idx1[0]])
        string += '%10.5f' % (varray[idx2[0]])
        string += '    1   -4     '
        string += '%10.8f' % (delv)
        string += '   12    1    1   13\n'
        fid.write(string)
        string = '%10.8f%10.5f' % (opd, varray[idx1[0]])
        string += '%10.5f' % (varray[idx2[0]])
        string += '    0   -4     '
        string += '%10.8f' % (delv)
        string += '   12    1    1   14\n'
        fid.write(string)
         
        #fid.write('%10.3f' % HWHM)
        #fid.write('%10.3f' % JV1)
        #fid.write('%10.3f' % JV2)
        #fid.write('%5i' % JEMIT)
        #fid.write('%5i' % JFNin)
        #fid.write('%5i' % MRATin)
        #fid.write('%10.3f' % DVOUT)
        #fid.write('%5i' % IUNIT)
        #fid.write('%5i' % IFILST)
        #fid.write('%5i' % NIFILS)
        #fid.write('%5i' % JUNIT)
        # Neglecting the IVX and NOFIX
        #fid.write('\n')
        fid.write('-1.\n')
        #end of Card 8.1

        blah = 0
        # This makes the TAPE27 and TAPE28
        # Settings taken from DDT rundecker.pro
        if ISCAN > 0 and IEMIT == 1:
            #Reinterpolate scanned data to regular coordinate system
            #Card 1.1
            fid.write('$ Transfer to ASCII plotting data (TAPES 27 and 28)\n')

            #Card 1.2
            fid.write(' HI=%1i' % 0)        #IHIRAC
            fid.write(' F4=%1i' % 0)        #ILBLF4
            fid.write(' CN=%1i' % 0)        #ICNTNM
            fid.write(' AE=%1i' % 0)        #IAERSL
            fid.write(' EM=%1i' % 0)        #IEMIT
            fid.write(' SC=%1i' % 0)        #ISCAN
            fid.write(' FI=%1i' % 0)        #IFILTR
            fid.write(' PL=%1i' % 1)        #IPLOT
            fid.write(' TS=%1i' % 0)        #ITEST
            fid.write(' AM=%1i' % 0)        #IATM
            fid.write(' MG=%1i' % 0)        #IMRG
            fid.write(' LA=%1i' % 0)        #ILAS
            fid.write(' MS=%1i' % 0)        #IOD
            fid.write(' XS=%1i' % 0)        #IXSECT
            fid.write('%5i' % 0)        #MPTS
            fid.write('%5i' % 0)        #NPTS
            fid.write('\n')
            fid.write("# Plot title not used\n")

            #Card 9.1 ..each card listed here will interpert TAPE13 and TAPE14 (from scanning)
            # Convert TAPE13
            string = '%10.5f%10.5f' % (varray[idx1[0]], varray[idx2[0]]) # V1 to V2
            string += '   10.2000  100.0000    5    0   13    0     1.000 0  0    0\n'
            string += '    0.0000    1.2000    7.0200    0.2000    4    0    1    1    0    0 0    3 27\n'
            # Convert TAPE14
            string += '%10.5f%10.5f' % (varray[idx1[0]], varray[idx2[0]]) # V1 to V2
            string += '   10.2000  100.0000    5    0   14    0     1.000 0  0    0\n'
            string += '    0.0000    1.2000    7.0200    0.2000    4    0    1    0    0    0 0    3 28\n'
            
            fid.write(string)
        
            """
            # Commented out to ensure the correct AERI stuff from rundecker.
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
            """
            fid.write('-1.\n')
            #end of Card 9.1  
            
    fid.write("%%%")
    print "Done."
