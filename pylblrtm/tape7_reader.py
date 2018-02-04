import numpy as np
import sys
"""
    SCRIPT NAME:
    tape7_reader.py

    SCRIPT AUTHOR:
    Greg Blumberg (OU/CIMMS)
    wblumberg@ou.edu

    DESCRIPTION:
    This script reads the output file named TAPE7 that is generated by the 
    LBLRTM.  The TAPE7 file contains information on the interpolated thermodynamic
    grid.  After reading in the file, the readTape function returns two variables.

    readTape accepts these arguments:
    filename - the filename/path to the TAPE7

    It returns these two arguments:
    height - the height grid
    temp - the temperature grid

    I honestly don't remember right now what the units are of these two...I guess
    I'll have to run it at some point.
"""

def readTape(filename):
    f = open(filename, 'r')
    line = f.readline()
    line = f.readline()

    i_form = long(line[0:2])
    n_layers = long(line[2:5])
    n_mols = long(line[8:10])
    secnt0 = float(line[10:20])
    h1 = float(line[40:48])
    h2 = float(line[52:60])
    ang = float(line[65:65+8])
    lens = long(line[78:80])
    print i_form, n_layers, n_mols
    print secnt0, h1, h2, ang, lens
    print line

    #Loop over layers 

    pressure = np.empty(n_layers)
    temperature = np.empty(n_layers)

    moldens = np.empty((n_mols, n_layers))
    #moldens_temp = np.empty(n_mols + 1 > 8)
    if n_mols <= 7:
        moldens_loc = np.arange(n_mols)
    else:
        moldens_loc = np.hstack((np.arange(7, np.int64), np.arange( n_mols - 7) + 8)) 

    path = np.arange(0, n_layers, 1, dtype=np.int64)
    p_top = np.arange(0, n_layers, 1, dtype=np.float64)
    p_bot = np.arange(0, n_layers, 1, dtype=np.float64)
    t_top = np.arange(0, n_layers, 1, dtype=np.float64)
    t_bot = np.arange(0, n_layers, 1, dtype=np.float64)
    z_top = np.arange(0, n_layers, 1, dtype=np.float64)
    z_bot = np.arange(0, n_layers, 1, dtype=np.float64)
    broaddens = np.arange(0, n_layers, 1, dtype=np.float64)
    
    ##Read and parse the first line as the format is a big different
    c_buf = f.readline()

    pressure[0] = float(c_buf[0:15])
    temperature[0] = float(c_buf[15:25])
    path[0] = int(c_buf[39:40])
    z_bot[0] - float(c_buf[41:48])
    p_bot[0] = float(c_buf[48:48+8])
    t_bot[0] = float(c_buf[56:56+7])
    z_top[0] = float(c_buf[63:63+7])
    p_top[0] = float(c_buf[70:78])
    t_top[0] = float(c_buf[78:78+7])

    moldens_tmp = np.asarray(f.readline().split(), dtype=np.float64)
    moldens[:,0] = moldens_tmp[moldens_loc]
    broaddens[0] = moldens_tmp[7]

    for i in np.arange(1, n_layers):
        c_buf = f.readline()
        pressure[i] = float(c_buf[0:15])
        temperature[i] = float(c_buf[15:25])
        path[i] = int(c_buf[39:40])
        z_top[i] = float(c_buf[63:70])
        t_top[i] = float(c_buf[78:78+7])
        p_top[i] = float(c_buf[70:78])
        p_bot[i] = p_top[i-1]
        t_bot[i] = t_top[i-1]
        z_bot[i] = z_top[i-1]

        moldens_tmp = np.asarray(f.readline().split())
        moldens[:,i] = moldens_tmp[moldens_loc]
        broaddens[i] = moldens_tmp[7]
    
    f.close()
    #print moldens
    return np.hstack((z_bot[0],z_top)), np.hstack((t_bot[0], t_top))
    
#readTape('../lblex/out/TAPE7')


