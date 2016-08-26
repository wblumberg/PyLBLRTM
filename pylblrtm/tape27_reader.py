import numpy as np

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
w1, rad1 = readTAPE27('/Users/greg.blumberg/clamps_calib/temp/TAPE27')
w2, rad2 = readTAPE27('/Users/greg.blumberg/clamps_calib/dave_lbl/TAPE27')
print np.mean(w1-w2)
print np.max(w1-w2)
print np.min(w1-w2)



print np.mean(rad1-rad2)
print np.max(rad1-rad2)
print np.min(rad1-rad2)


