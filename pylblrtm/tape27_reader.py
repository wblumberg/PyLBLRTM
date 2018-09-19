import numpy as np
from pylab import *
import sys

def readTAPE27(filename, scale=1e7):
    fn = open(filename, 'r')
    rad = []
    wnum = []

    for i, line in enumerate(fn.readlines()):
        if i > 26:
            parsed = line.split()
            rad.append(parsed[1])
            wnum.append(parsed[0])
    fn.close()
    return np.asarray(wnum, dtype=float), np.asarray(rad, dtype=float)*scale
w1, rad1 = readTAPE27(sys.argv[1])
plot(w1, rad1)
xlabel("Wavenumber [cm-1]")
ylabel("Radiance [RU]")
show()


