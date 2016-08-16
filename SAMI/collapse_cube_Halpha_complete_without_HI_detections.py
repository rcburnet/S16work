#Script to create collapsed cube between FWHM of Halpha profile (ie. create Halpha image) and plot spectrum of 58 SAMI targets that don't have HI detections. Note that the flux is in units of 10**(-16) erg /s /cm**2 /angstrom /pixel, not in proper units. May need to multiply by number of pixels to get right units (?)

from astropy.io import fits
import numpy as np
import os
import fnmatch
import matplotlib.pyplot as plt


hdulist = []
filename_list = []

for root,dirnames,filenames in os.walk('/home/rburnet/SAMI/data/without_HI_detections/'):
    for filename in fnmatch.filter(filenames, '[0-9]*red*[0-9].fits'):
        hdulist.append(os.path.join(root, filename))
        filename_list.append(filename)

#Extract GAMA id's from filenames
GAMA_name_list_no_HI_detections = []

for i in range(len(filename_list)):
    if '_' in filename_list[i][0:6]:
        GAMA_name_list_no_HI_detections.append(filename_list[i][0:5])
    else:
        GAMA_name_list_no_HI_detections.append(filename_list[i][0:6])


#Extract required data (redshifts)
sami = open('/home/rburnet/S16work/SAMI/SAMI_EarlyDataRelease_modified.txt')
sami_lines = sami.readlines()

for i in range(len(sami_lines)):
    sami_lines[i] = sami_lines[i].split(' ')

sami.close()

sami_z = []

for j in range(len(GAMA_name_list_no_HI_detections)):
    for i in range(len(sami_lines)):
        if i != 0 and i != 108:
            if GAMA_name_list_no_HI_detections[j] in sami_lines[i]:
                sami_z.append(float(sami_lines[i][6]))

for i in range(len(hdulist)):
    hdulist1 = fits.open(hdulist[i])

    data = hdulist1[0].data

    data = np.nan_to_num(data)

    for j in range(len(GAMA_name_list_no_HI_detections)):
        if GAMA_name_list_no_HI_detections[j] in hdulist[i]:
            z = sami_z[j]

    centre_wavelength = hdulist1[0].header['CRVAL3']
    increment = hdulist1[0].header['CDELT3']
    x = np.linspace(centre_wavelength - increment * 1023, centre_wavelength + increment * 1024, 2048)

    halpha = 6562.8 #Angstroms

    l = halpha*(1+z)
    NII_1 = 6548 * ( 1 + z )
    NII_2 = 6584 * ( 1 + z )

    Halpha_flux = []    #the flux within the Halpha profile
    continuum_flux = [] #the continuum flux (between Halpha+7 and NII_2-7)

    # retrieve continuum flux
    data = hdulist1[0].data
    for j in range(len(data)):
        if x[j] > l+7 and x[j] < NII_2-7:    #Only place values between Halpha+7 and NII_2-7 (to get continuum)
            continuum_flux.append(data[j])

    #continuum is mean over continuum wavelength range for each pixel, 50x50 array where each pixel is the continuum flux mean for that specific pixel.
    continuum_flux = np.mean(continuum_flux, axis=0)

    # subtract Halpha profile from continuum, then integrate over Halpha profile to get Halpha flux
    for j in range(len(data)):
        if x[j] > l-7 and x[j] < l+7:   #look over Halpha profile and a little in the continuum to subtract continuum from flux and then integrate
            Halpha_flux.append(data[j] - continuum_flux)

    Halpha_flux = np.array(Halpha_flux)
    
    summed_flux = np.sum(Halpha_flux, axis=0)

    hdulist1[0].data = summed_flux

    hdulist1.writeto(filename_list[i]+'collapsed.fits')

    hdulist1.close()
