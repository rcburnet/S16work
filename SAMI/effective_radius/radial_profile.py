# Script to plot radial profile of fits file.

from astropy.io import fits
import os
import fnmatch
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator

minorLocator = AutoMinorLocator()

def radial_profile(data, center):
    x, y = np.indices((data.shape))
    r = np.sqrt((x - center[0])**2 + (y - center[1])**2)
    r = r.astype(np.int)

    tbin = np.bincount(r.ravel(), data.ravel())
    nr = np.bincount(r.ravel())
    radialprofile = tbin / nr
    return radialprofile 

hdulist = []
filename_list = []

for root,dirnames,filenames in os.walk('/home/rburnet/SAMI/data/Halpha/without_HI_detections/'): # can change to without_HI_detections
    for filename in fnmatch.filter(filenames, '*.fits'):
        hdulist.append(os.path.join(root, filename))
        filename_list.append(filename)

#Extract GAMA id's from filenames
GAMA_name_list = []

for i in range(len(filename_list)):
    if '_' in filename_list[i][0:6]:
        GAMA_name_list.append(filename_list[i][0:5])
    else:
        GAMA_name_list.append(filename_list[i][0:6])

for k in range(len(hdulist)):
    
    print filename_list[k]
    fitsFile = fits.open(hdulist[k])
    img = fitsFile[0].data
    img[np.isnan(img)] = 0

    #center = np.unravel_index(img.argmax(), img.shape)
    center = (fitsFile[0].header['CRPIX2'], fitsFile[0].header['CRPIX1'])
    rad_profile = radial_profile(img, center)
    I_0 = rad_profile[0]

    I_e = I_0 * np.e**(-1.678) #for spiral galaxies

    fig, ax = plt.subplots()
    plt.plot(rad_profile[0:22], 'x-')

    I_e_list = []
    for i in range(25):
        I_e_list.append(I_e)

    x = np.linspace(0,25,25)

    plt.plot(x, I_e_list)    #plot I_e

    #attempt at finding point of intersection
    '''
    i = 1
    j = 0
    while i>0:
        i = rad_profile[j] - I_e
        j += 1

    print j-1   #R at point just past intersection point

    print rad_profile[j-1]  #I at point just past intersection point

    slope = (rad_profile[j-2] - rad_profile[j-1])   #slope of line between the points on either side of the intersection point
    '''
    ax.xaxis.set_minor_locator(minorLocator)

    plt.tick_params(which='both', width=2)
    plt.tick_params(which='major', length=7)
    plt.tick_params(which='minor', length=4, color='r')
    plt.grid()
    #ax.set_ylabel(fitsFile[0].header['Label'] + " (" + fitsFile[0].header['BUNIT'] + ")")
    ax.set_ylabel("Flux (1e-16 Ergs/s/cm$^2$)")
    ax.set_xlabel("Radius (Pixels)")
    plt.title("Radial Profile of "+str(GAMA_name_list[k]))
    plt.grid(which="minor")
    plt.savefig('figures/'+str(GAMA_name_list[k])+'_radial_profile.pdf')
    plt.close()

    fitsFile.close()
