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

### Read data cubes
hdulist_final = []
hdulist_cubes = []
image_list = []

for root,dirnames,filenames in os.walk('/home/rburnet/reflex/project_data/OBs_combined/'):
    for filename in fnmatch.filter(filenames, '*collapsed*.fits'):
        hdulist_final.append(os.path.join(root, filename))

for root,dirnames,filenames in os.walk('/home/rburnet/reflex/project_data/OBs_combined/'):
    for filename in fnmatch.filter(filenames, '*combine_sci_reconstructed*.fits'):
        hdulist_cubes.append(os.path.join(root, filename))

for root,dirnames,filenames in os.walk('/home/rburnet/S16work/diagnostic_of_objects/figures/HST_images/'):
    for filename in fnmatch.filter(filenames, 'CL003[0-9]-Target_*.png'):
        image_list.append(os.path.join(root, filename))

hdulist_final = sorted(hdulist_final)
hdulist_cubes = sorted(hdulist_cubes)
image_list = sorted(image_list)

### Read txt table files for target/arm info
file_CL0034 = open('../txt_tables/KMOS_GCLASS_CLUS0034.txt','r')
file_CL0036_YJ = open('../txt_tables/KMOS_GCLASS_CLUS0036_YJ.txt','r')
file_CL0036_IZ = open('../txt_tables/KMOS_GCLASS_CLUS0036_IZ.txt','r')

CL0034_info = file_CL0034.readlines()
CL0036_YJ_info = file_CL0036_YJ.readlines()
CL0036_IZ_info = file_CL0036_IZ.readlines()

file_CL0034.close()
file_CL0036_YJ.close()
file_CL0036_IZ.close()

for i in range(len(CL0034_info)):
    CL0034_info[i] = CL0034_info[i].split(' ')
    for j in range(len(CL0034_info[i])):
        if j < 3:
            CL0034_info[i][j] = int(CL0034_info[i][j])
        else:
            CL0034_info[i][j] = float(CL0034_info[i][j])

for i in range(len(CL0036_YJ_info)):
    CL0036_YJ_info[i] = CL0036_YJ_info[i].split(' ')
    for j in range(len(CL0036_YJ_info[i])):
        if j < 3:
            CL0036_YJ_info[i][j] = int(CL0036_YJ_info[i][j])
        else:
            CL0036_YJ_info[i][j] = float(CL0036_YJ_info[i][j])

for i in range(len(CL0036_IZ_info)):
    CL0036_IZ_info[i] = CL0036_IZ_info[i].split(' ')
    for j in range(len(CL0036_IZ_info[i])):
        if j < 3:
            CL0036_IZ_info[i][j] = int(CL0036_IZ_info[i][j])
        else:
            CL0036_IZ_info[i][j] = float(CL0036_IZ_info[i][j])

x = np.linspace(0.78, 1.09, 2048)

### Now plot stuff

for i in range(len(hdulist_final)):
    hdufits_final = fits.open(hdulist_final[i])
    hdufits_cubes = fits.open(hdulist_cubes[i])
    cube_data = hdufits_cubes[1].data
    img_data = hdufits_final[1].data
    hdufits_final.close()
    hdufits_cubes.close()

    start = len('/home/rburnet/reflex/project_data/OBs_combined/CL0034_YJ/')
    
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
                if 'CL0034' in image_list[j]:
                    if '_' in image_list[j][-6:-4]:
                        target_name = int(image_list[j][-5:-4])
                    else:
                        target_name = int(image_list[j][-6:-4])
                    #print image_list[j][-7:-5]
                    for k in range(len(CL0034_info)):
                        if target_name == CL0034_info[k][0]:
                            mag_exp = CL0034_info[k][7]

                if 'CL0036' in image_list[j]:
                    if '_' in image_list[j][-6:-4]:
                        target_name = int(image_list[j][-5:-4])
                    else:
                        target_name = int(image_list[j][-6:-4])
                    for k in range(len(CL0036_YJ_info)):    #doesn't matter if you use YJ or IZ since target names are the same between them. Only the arms names are different.
                        if target_name == CL0036_YJ_info[k][0]:
                            mag_exp = CL0036_YJ_info[k][7]
                #else:
                #    mag_exp = 0.00
                plt.text(18,10,'m$_{expected}$ = '+'%.2f' % round(mag_exp,2))

                plt.savefig('figures/'+hdulist_final[i][-35:]+'.png', bbox_inches = 'tight')
                plt.clf()
                plt.close()

                im = Image.open('figures/'+hdulist_final[i][-35:]+'.png')
                im.crop([40,105,800,410]).save('figures/'+hdulist_final[i][-35:]+'.png', 'PNG')
                im.close()
