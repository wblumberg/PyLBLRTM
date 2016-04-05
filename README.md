##PyLBLRTM

PyLBLRTM is a collection of Python scripts to ease the use of the line-by-line radiative transfer model (LBLRTM) managed by AER.  These scripts allow for the following:

1.) Creation of the input file for the LBLRTM (TAPE5).  Different configurations of the TAPE5 exist, but the writer contained in this package allows for both the creation of optical depth files (needed for LBLDIS and any of your own radiative transfer calculations) and radiance calculations in the TAPE12 file.

2.) Reading of the TAPE7 file that contains information about the grid the LBLRTM is using to compute optical depths, etc.

3.) The ability to read in "panel" files, which are binary formatted files output by the LBLRTM.

4.) The ability to read in a directory of files output from the LBLRTM script "lblrun".  The data from the LBLRTM are contained within an object called LBLPkg.  From LBLPkg, upwelling and downwelling calculations of radiance can be performed.  LBLPkg does not support including scattering in the radiative transfer calculations.  If you want scattering, this isn't the package for it, bub.

-

These files were used in October 2013 to run the LBLRTM and analyze the output for AERI files.
They were written by Greg Blumberg (OU/CIMMS).



