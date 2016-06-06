#!/usr/bin/python

####
# script to present collapsed cubes of object, image of object from HST images, and calculate/expected magnitudes into a visual format for diagnostic, to see if the objects in the raw data agree with the HST images. To produce final diagnostic report, go to /report/ and run diagnostic_report_OB(1|2).tex
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
from PIL import Image

hdulist_final = []
hdulist_cubes = []
image_list = []

for root,dirnames,filenames in os.walk('/home/rburnet/reflex/project_data/reflex_end_products/2016-05-24T13:39:06/'):
    for filename in fnmatch.filter(filenames, 'CL003[0-9]-IZ-OB-[1|2]_COMBINE_SCI_RECONSTRUCTED_COLL_*.fits'):
        hdulist_final.append(os.path.join(root, filename))

for root,dirnames,filenames in os.walk('/home/rburnet/reflex/project_data/reflex_end_products/2016-05-24T13:39:06/'):
    for filename in fnmatch.filter(filenames, 'CL003[0-9]-IZ-OB-[1|2]_COMBINE_SCI_RECONSTRUCTED_[0-9]*.fits'):
        hdulist_cubes.append(os.path.join(root, filename))

for root,dirnames,filenames in os.walk('/home/rburnet/S16work/diagnostic_of_objects/figures/HST_images/'):
    for filename in fnmatch.filter(filenames, 'CL003[0-9]-Target_*.png'):
        image_list.append(os.path.join(root, filename))

hdulist_final = sorted(hdulist_final)
hdulist_cubes = sorted(hdulist_cubes)
image_list = sorted(image_list)

x = np.linspace(0.78, 1.09, 2048)

for i in range(len(hdulist_final)):
    hdufits_final = fits.open(hdulist_final[i])
    hdufits_cubes = fits.open(hdulist_cubes[i])
    cube_data = hdufits_cubes[1].data
    img_data = hdufits_final[1].data
    hdufits_final.close()
    hdufits_cubes.close()

    start = len('/home/rburnet/reflex/project_data/reflex_end_products/2016-05-24T13:39:06/KMOS.2014-07-11T06:34:52.380_tpl/')
    
    for j in range(len(image_list)):
        if hdulist_final[i][start:start+6] == image_list[j][63:69]:
            if hdulist_final[i][-7:-5] == image_list[j][-6:-4]:
                width = img_data.shape[0]
                height = img_data.shape[1]

                img_data = np.array(img_data, dtype=float)
                img_data = np.nan_to_num(img_data)
                min_val = np.min(img_data)

                index = (792, 793) #index of 900nm slice
                array1 = np.nan_to_num(cube_data[792])
                array2 = np.nan_to_num(cube_data[793])
                array1 = np.sum(array1)
                array2 = np.sum(array2)
                flux = np.mean([array1, array2])
                mag = -2.5*np.log10(3.34e4*(9000.0**2)*flux)+8.90

                fig = plt.figure()
                plt.axis('off')
                a = fig.add_subplot(121)
                a.axis('off')
                new_img = img_scale.linear(img_data, scale_min=min_val)
                imgplot = plt.imshow(new_img, interpolation='nearest', origin='lower', cmap='gray')
                if flux > 0:
                    plt.text(18,13,'m$_{calculated}$ = '+'%.2f' % round(mag,2))
                else:
                    plt.text(18,13,'m$_{calculated}$ = Not computed, flux negative')
                if 'CL0034-Target_5.' in image_list[j]:
                    mag_exp = 20.46
                if 'CL0034-Target_6.' in image_list[j]:
                    mag_exp = 22.37
                if 'CL0034-Target_7.' in image_list[j]:
                    mag_exp = 21.06
                if 'CL0034-Target_11' in image_list[j]:
                    mag_exp = 19.43
                if 'CL0034-Target_12' in image_list[j]:
                    mag_exp = 20.55
                if 'CL0034-Target_14' in image_list[j]:
                    mag_exp = 22.49
                if 'CL0034-Target_16' in image_list[j]:
                    mag_exp = 22.50
                if 'CL0034-Target_20' in image_list[j]:
                    mag_exp = 21.73
                if 'CL0034-Target_21' in image_list[j]:
                    mag_exp = 22.91
                if 'CL0034-Target_22' in image_list[j]:
                    mag_exp = 22.27
                if 'CL0034-Target_23' in image_list[j]:
                    mag_exp = 22.08
                if 'CL0034-Target_26' in image_list[j]:
                    mag_exp = 22.40
                if 'CL0034-Target_30' in image_list[j]:
                    mag_exp = 20.03
                if 'CL0034-Target_31' in image_list[j]:
                    mag_exp = 20.84
                if 'CL0036-Target_3.' in image_list[j]:
                    mag_exp = 20.78
                if 'CL0036-Target_4.' in image_list[j]:
                    mag_exp = 20.24
                if 'CL0036-Target_5.' in image_list[j]:
                    mag_exp = 20.80
                if 'CL0036-Target_6.' in image_list[j]:
                    mag_exp = 21.89
                if 'CL0036-Target_7.' in image_list[j]:
                    mag_exp = 21.53
                if 'CL0036-Target_8.' in image_list[j]:
                    mag_exp = 21.15
                if 'CL0036-Target_10' in image_list[j]:
                    mag_exp = 22.06
                if 'CL0036-Target_17' in image_list[j]:
                    mag_exp = 21.67
                if 'CL0036-Target_19' in image_list[j]:
                    mag_exp = 21.68
                if 'CL0036-Target_26' in image_list[j]:
                    mag_exp = 19.27
                if 'CL0036-Target_28' in image_list[j]:
                    mag_exp = 20.35
                if 'CL0036-Target_30' in image_list[j]:
                    mag_exp = 21.30
                if 'CL0036-Target_34' in image_list[j]:
                    mag_exp = 20.36
                #else:
                #    mag_exp = 0.00
                plt.text(18,10,'m$_{expected}$ = '+'%.2f' % round(mag_exp,2))

                plt.savefig('figures/'+hdulist_final[i][-53:]+'.png', bbox_inches = 'tight')
                plt.clf()
                plt.close()

                im = Image.open('figures/'+hdulist_final[i][-53:]+'.png')
                im.crop([40,105,800,410]).save('figures/'+hdulist_final[i][-53:]+'.png', 'PNG')
                im.close()
