#!/usr/bin/python

####
#script to modify obs mode in header of raw fits files
#for some reason, script won't run and write fits files unless you copy and past this code into python terminal. I don't understand why. It's really annoying.
####

import numpy as np
import scipy as sp
import glob
import os
import fnmatch
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy.io import fits

hdulist=[]

for root,dirnames,filenames in os.walk('/home/rburnet/reflex/project_data/reflex_input/kmos/data_with_raw_calibs/'):
    for filename in fnmatch.filter(filenames, '*.fits'):
        hdulist.append(os.path.join(root, filename))

for i in range(len(hdulist)):
    hdufits = fits.open(hdulist[i])
    hdufits[0].header
    try: hdufits[0].header['HIERARCH ESO TPL MODE OBS'];hdufits[0].header['HIERARCH ESO TPL MODE OBS'] = 'STARE';print 'worked'
    except: print hdulist[i]+' has no relevant header for NOD TO SKY'
    hdufits.writeto(hdulist[i][0:52]+'STARE/'+hdulist[i][73:])
    hdufits.close()
