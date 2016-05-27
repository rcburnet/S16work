#!/usr/bin/python

####
# script to calculate f = <T_lambda(sky)>/<S_lambda(sky)> , the fraction of average flux of the target frame sky exposures to the average flux of the sky frame sky exposures to see if they agree.
####

import numpy as np
import scipy as sp
import glob
import os
import fnmatch
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy.io import fits

hdulist_sci = []
hdulist_sky = []

for root,dirnames,filenames in os.walk('/home/rburnet/reflex/project_data/reflex_end_products/'):
    for filename in fnmatch.filter(filenames, 'CL003[0-9]-YJ-OB-[1|2]_COMBINE_SCI_RECONSTRUCTED_ARM*_SCI.fits'):
        hdulist_sci.append(os.path.join(root, filename))
        hdulist_sky.append(os.path.join(root, filename[:-8]+'SKY.fits'))

f = []  #inital fraction for each wavelength slice of every pair of files. Has shape (#of pairs of files, 2048)

for i in range(len(hdulist_sci)):
    f.append([])
    hdufits_sci = fits.open(hdulist_sci[i])
    hdufits_sky = fits.open(hdulist_sky[i])
    for j in range(len(hdufits_sci[1].data)):
        avg_flux_sci = np.mean(hdufits_sci[1].data[j][~np.isnan(hdufits_sci[1].data[j])])
        avg_flux_sky = np.mean(hdufits_sky[1].data[j][~np.isnan(hdufits_sky[1].data[j])])
        if abs(avg_flux_sci/avg_flux_sky) == np.inf:    #this can occur if there is a slice with data in the sci file, but no data in the corresponding slice in the sky file
            f[i].append(np.nan)                         #if this occurs, append NaN (to be ignored later on)
        else:
            f[i].append(avg_flux_sci/avg_flux_sky)
    hdufits_sci.close()
    hdufits_sky.close()

f2 = [] #f with NaN values removed
f3 = [] #average f for every file with NaN removed
f = np.array(f)

for i in range(len(f)):
    f2_i = f[i]
    f2.append(f2_i[~np.isnan(f2_i)])
    f3.append(np.mean(f2[i]))


print 'f from all files = ',f3
print 'average f for all files = ',np.mean(f3)

x = np.linspace(1,1.35882,2048)

#loop to turn all NaN values in f array to 1.0 (expected value)
for i in range(len(f)):
    for j in range(len(f[i])):
        if f[i][j] != f[i][j]:
            f[i][j] = 1.0
    #now plot the resulting f's for every file and every wavelength slice
    if not os.path.exists('./figures/laptop/target_sky_fraction/'):
        os.makedirs('./figures/laptop/target_sky_fraction/')
    line1, = plt.plot(x,f[i])
    line5, = plt.plot((1.0,1.35882),(1.0,1.0))
    plt.title('sci_flux/sky_flux from sky arms vs Wavelength from \n '+hdulist_sci[i]+' \n and \n  '+hdulist_sky[i],fontsize=10)
    plt.xlabel('Wavelength (microns)',fontsize=10)
    plt.ylabel('sci_flux/sky_flux',fontsize=10)
    plt.axis('tight')
    plt.ylim([0,2])
    fig = plt.gcf()
    fig.set_size_inches(18.5, 5.5)
    plt.setp(line1, linewidth=0.1, color='b')
    plt.setp(line5, linewidth=0.1)
    print './figures/laptop/target_sky_fraction/'+hdulist_sci[i][-55:]+'.pdf'
    plt.savefig('./figures/laptop/target_sky_fraction/'+hdulist_sci[i][-55:]+'.pdf')
    plt.close()
