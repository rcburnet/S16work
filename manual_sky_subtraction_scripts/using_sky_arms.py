#!/usr/bin/python

import numpy as np
import scipy as sp
import glob
import os
import fnmatch
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy.io import fits

hdulist_Ssci = []
hdulist_sky = []
hdulist = []

for root,dirnames,filenames in os.walk('/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/'):
    for filename in fnmatch.filter(filenames, 'CL003[4|6]-[Y|I][J|Z]-OB-[1|2]_SCI_RECONSTRUCTED_KMOS*'):
        hdulist.append(os.path.join(root, filename))

hdulist = sorted(hdulist)

fitsfiles = {}
for i in range(len(hdulist)):
    fitsfiles[i] = fits.open(hdulist[i])    #fitsfiles is a dictionary, each dict key has a value which is an open fits file. Each fits file is 25 elements long (1 for each IFU, 1 primary header)

#0 - 4 is CL0034 YJ OB1 ABAAB, 5 - 7 is CL0034 YJ OB2 ABA, 8 - 12 is CL0036 YJ OB1 ABAAB, 13 - 15 is CL0036 YJ OB2 ABA, 16 - 20 is CL034 IZ OB1 ABAAB, 21 - 25 is CL0034 IZ OB2 ABAAB, 26 - 30 is CL0036 IZ OB1 ABAAB, 31 - 35 is CL0036 IZ OB2 ABAAB

Sky_Arm_2 = np.array(fitsfiles[1][2].data)
Sci_Arm_2 = np.array(fitsfiles[0][2].data)
Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
Sci_Arm_sub = []
for i in range(len(Sky_Arm_2)):
    Sci_Arm_sub.append([])
    Sky_flux = np.sum(Sky_Arm_2[i])
    Sci_flux = np.sum(Sci_Arm_2[i])
    f = Sci_flux / Sky_flux
    new_flux = Sci_Arm_2[i] - f * Sky_Arm_2[i]
    Sci_Arm_sub[i]  = new_flux

Sci_Arm_sub = np.array(Sci_Arm_sub)
Sci_Arm_sub = np.nan_to_num(Sci_Arm_sub)
print np.sum(Sci_Arm_sub)

'''
CL0034_YJ_OB1_fits_A1 = fits.open(hdulist[0])
CL0034_YJ_OB1_fits_B1 = fits.open(hdulist[1])
CL0034_YJ_OB1_fits_A2 = fits.open(hdulist[2])
CL0034_YJ_OB1_fits_A3 = fits.open(hdulist[3])
CL0034_YJ_OB1_fits_B2 = fits.open(hdulist[4])
CL0034_YJ_OB2_fits_A1 = fits.open(hdulist[5])
CL0034_YJ_OB2_fits_B1 = fits.open(hdulist[6])
CL0034_YJ_OB2_fits_A2 = fits.open(hdulist[7])
CL0036_YJ_OB1_fits_A1 = fits.open(hdulist[8])
CL0036_YJ_OB1_fits_B1 = fits.open(hdulist[9])
CL0036_YJ_OB1_fits_A2 = fits.open(hdulist[10])
CL0036_YJ_OB1_fits_A3 = fits.open(hdulist[11])
CL0036_YJ_OB1_fits_B2 = fits.open(hdulist[12])
CL0036_YJ_OB2_fits_A1 = fits.open(hdulist[13])
CL0036_YJ_OB2_fits_B1 = fits.open(hdulist[14])
CL0036_YJ_OB2_fits_A2 = fits.open(hdulist[15])
CL0034_IZ_OB1_fits_A1 = fits.open(hdulist[16])
CL0034_IZ_OB1_fits_B1 = fits.open(hdulist[17])
CL0034_IZ_OB1_fits_A2 = fits.open(hdulist[18])
CL0034_IZ_OB1_fits_A3 = fits.open(hdulist[19])
CL0034_IZ_OB1_fits_B2 = fits.open(hdulist[20])
CL0034_IZ_OB2_fits_A1 = fits.open(hdulist[21])
CL0034_IZ_OB2_fits_B1 = fits.open(hdulist[22])
CL0034_IZ_OB2_fits_A2 = fits.open(hdulist[23])
CL0034_IZ_OB2_fits_A3 = fits.open(hdulist[24])
CL0034_IZ_OB2_fits_B2 = fits.open(hdulist[25])
CL0036_IZ_OB1_fits_A1 = fits.open(hdulist[26])
CL0036_IZ_OB1_fits_B1 = fits.open(hdulist[27])
CL0036_IZ_OB1_fits_A2 = fits.open(hdulist[28])
CL0036_IZ_OB1_fits_A3 = fits.open(hdulist[29])
CL0036_IZ_OB1_fits_B2 = fits.open(hdulist[30])
CL0036_IZ_OB2_fits_A1 = fits.open(hdulist[31])
CL0036_IZ_OB2_fits_B1 = fits.open(hdulist[32])
CL0036_IZ_OB2_fits_A2 = fits.open(hdulist[33])
CL0036_IZ_OB2_fits_A3 = fits.open(hdulist[34])
CL0036_IZ_OB2_fits_B2 = fits.open(hdulist[35])
'''
