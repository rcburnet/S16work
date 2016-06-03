#!/usr/bin/python

####
# script to present collapsed cubes of object, image of object from HST images, and calculate/expected magnitudes into a visual format for diagnostic, to see if the objects in the raw data agree with the HST images.
####

import numpy as np
import scipy as sp
import glob
import os
import fnmatch
import matplotlib.pyplot as plt
import matplotlib.image as mpimg
from matplotlib import rc
from astropy.io import fits
import img_scale
import pylab

hdulist_final = []

for root,dirnames,filenames in os.walk('/home/rburnet/reflex/project_data/reflex_end_products/2016-05-24T13:39:06/'):
    for filename in fnmatch.filter(filenames, 'CL003[0-9]-[Y|I][J|Z]-OB-[1|2]_COMBINE_SCI_RECONSTRUCTED_COLL_*.fits'):
        hdulist_final.append(os.path.join(root, filename))


for i in range(len(hdulist_final)):
    hdufits_final = fits.open(hdulist_final[i])
    img_data = hdufits_final[1].data
    hdufits_final.close()

    width = img_data.shape[0]
    height = img_data.shape[1]

    img_data = np.array(img_data, dtype=float)
    img_data = np.nan_to_num(img_data)
    min_val = np.min(img_data)

    fig = plt.figure()
    a = fig.add_subplot(1,2,1)
    new_img = img_scale.linear(img_data, scale_min=min_val)
    imgplot = plt.imshow(new_img, interpolation='nearest', origin='lower', cmap='gray')
    a.set_title('Linear')
    plt.colorbar(orientation='horizontal')

    a = fig.add_subplot(1,2,2)
    new_img = img_scale.sqrt(img_data, scale_min=min_val)
    imgplot = plt.imshow(new_img, interpolation='nearest', origin='lower', cmap='gray')
    a.set_title('Sqrt')
    plt.colorbar(orientation='horizontal')

    plt.savefig('figures/'+hdulist_final[i][-53:]+'.png')
    plt.clf()
    plt.close()
