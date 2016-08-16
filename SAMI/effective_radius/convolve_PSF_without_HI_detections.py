# Convolve PSF and image fits of 58 SAMI targets without HI detections

from astropy.convolution import convolve, convolve_fft
from astropy.io import fits
import os
import fnmatch
import numpy as np

hdulist = []
filename_list = []
secondary_star_list = []

for root,dirnames,filenames in os.walk('/home/rburnet/SAMI/data/Halpha/without_HI_detections/'):
    for filename in fnmatch.filter(filenames, '*.fits'):
        hdulist.append(os.path.join(root, filename))
        filename_list.append(filename)

for root,dirnames,filenames in os.walk('/home/rburnet/SAMI/data/standard_stars/'):
    for filename in fnmatch.filter(filenames, '1*.fits'):
        secondary_star_list.append(os.path.join(root, filename))

for i in range(len(hdulist)):

    image = fits.open(hdulist[i])

    std_name = image[0].header['STDNAME']
    plate_id = image[0].header['PLATEID']
    print filename_list[i]
    
    for j in range(len(secondary_star_list)):
        if std_name in secondary_star_list[j]:
            psf_image = fits.open(secondary_star_list[j])

    image_data = image[0].data
    psf_image_data = psf_image[0].data[100:2000]
    psf_image_data = np.nanmean(psf_image_data, axis=0)

    convolved_data = convolve_fft(image_data, psf_image_data)

    for j in range(len(image_data)):
        for k in range(len(image_data[j])):
            if image_data[j][k] != image_data[j][k]:
                convolved_data[j][k] = np.nan

    image[0].data = convolved_data

    image.writeto(filename_list[i][:-19]+'_PSF_convolved.fits')

    image.close()
    psf_image.close()
