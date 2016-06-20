#!/usr/bin/python

####
# Script that will use sky arms from science frame and it's corresponding sky arm from sky frame to properly scale flux and then carry out sky subtraction
####

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

for i in range(len(fitsfiles)):
    #CL0034 AB pairs
    if i == 0 or i == 3 or i == 5 or i == 16 or i == 19 or i == 21 or i == 24:
        Sky_Arm_2 = np.array(fitsfiles[i+1][2].data)
        Sci_Arm_2 = np.array(fitsfiles[i][2].data)
        Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
        Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
        f1_list = []
        for k in range(len(Sky_Arm_2)):
            f1_list.append([])
            Sky_flux = np.mean(Sky_Arm_2[k])
            Sci_flux = np.mean(Sci_Arm_2[k])
            f = Sci_flux / Sky_flux
            f1_list[k]  = f

        if i < 15:
            Sky_Arm_2 = np.array(fitsfiles[i+1][5].data)
            Sci_Arm_2 = np.array(fitsfiles[i][5].data)
            Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
            Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
            f2_list = []
            for k in range(len(Sky_Arm_2)):
                f2_list.append([])
                Sky_flux = np.mean(Sky_Arm_2[k])
                Sci_flux = np.mean(Sci_Arm_2[k])
                f = Sci_flux / Sky_flux
                f2_list[k]  = f
        else:
            f2_list = f1_list

        for k in range(len(fitsfiles[i])):
            try:
                Sky_Arm_2 = np.array(fitsfiles[i+1][k].data)
                Sci_Arm_2 = np.array(fitsfiles[i][k].data)
                Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
                Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
                Sci_Arm_sub1 = []
                Sci_Arm_sub2 = []
                for j in range(len(Sky_Arm_2)):
                    Sci_Arm_sub1.append([])
                    Sci_Arm_sub2.append([])
                    new_flux1 = Sci_Arm_2[j] - f1_list[j] * Sky_Arm_2[j]
                    new_flux2 = Sci_Arm_2[j] - f2_list[j] * Sky_Arm_2[j]
                    Sci_Arm_sub1[j] = new_flux1
                    Sci_Arm_sub2[j] = new_flux2
                sub1 = np.mean(Sci_Arm_sub1)
                sub2 = np.mean(Sci_Arm_sub2)
                if sub1 > sub2:
                    Sci_Arm_sub = Sci_Arm_sub1
                else:
                    Sci_Arm_sub = Sci_Arm_sub2
                Sci_Arm_sub = np.array(Sci_Arm_sub)
                fitsfiles[i][k].data = Sci_Arm_sub
            except:
                print "no data for arm "+str(k)
            #hdu = fits.PrimaryHDU(header=fitsfiles[i][0].header)
            #hdulist_file = fits.HDUList([hdu])
            #hdulist_file.append(fitsfiles[i][k])
            #hdulist_file.writeto('/home/rburnet/reflex/project_data/sky_sub/'+hdulist[i][-66:-5]+'_arm'+str(k)+'.fits')
        fitsfiles[i].writeto('/home/rburnet/reflex/project_data/sky_sub/'+hdulist[i][-66:])

    #CL0036 AB pairs
    if i == 8 or i == 11 or i == 13 or i == 26 or i == 29 or i == 31 or i == 34:
        if i < 15:
            Sky_Arm_2 = np.array(fitsfiles[i+1][5].data)
            Sci_Arm_2 = np.array(fitsfiles[i][5].data)
        if i > 15:
            Sky_Arm_2 = np.array(fitsfiles[i+1][9].data)
            Sci_Arm_2 = np.array(fitsfiles[i][9].data)
        Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
        Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
        f1_list = []
        for k in range(len(Sky_Arm_2)):
            f1_list.append([])
            Sky_flux = np.mean(Sky_Arm_2[k])
            Sci_flux = np.mean(Sci_Arm_2[k])
            f = Sci_flux / Sky_flux
            f1_list[k]  = f

        if i < 15:
            Sky_Arm_2 = np.array(fitsfiles[i+1][17].data)
            Sci_Arm_2 = np.array(fitsfiles[i][17].data)
        if i > 15:
            Sky_Arm_2 = np.array(fitsfiles[i+1][13].data)
            Sci_Arm_2 = np.array(fitsfiles[i][13].data)
        Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
        Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
        f2_list = []
        for k in range(len(Sky_Arm_2)):
            f2_list.append([])
            Sky_flux = np.mean(Sky_Arm_2[k])
            Sci_flux = np.mean(Sci_Arm_2[k])
            f = Sci_flux / Sky_flux
            f2_list[k]  = f

        Sky_Arm_2 = np.array(fitsfiles[i+1][21].data)
        Sci_Arm_2 = np.array(fitsfiles[i][21].data)
        Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
        Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
        f3_list = []
        for k in range(len(Sky_Arm_2)):
            f3_list.append([])
            Sky_flux = np.mean(Sky_Arm_2[k])
            Sci_flux = np.mean(Sci_Arm_2[k])
            f = Sci_flux / Sky_flux
            f3_list[k]  = f


        if i < 15:
            Sky_Arm_2 = np.array(fitsfiles[i+1][13].data)
            Sci_Arm_2 = np.array(fitsfiles[i][13].data)
            Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
            Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
            f4_list = []
            for k in range(len(Sky_Arm_2)):
                f4_list.append([])
                Sky_flux = np.mean(Sky_Arm_2[k])
                Sci_flux = np.mean(Sci_Arm_2[k])
                f = Sci_flux / Sky_flux
                f4_list[k]  = f
        else:
            f4_list = f1_list


        for k in range(len(fitsfiles[i])):
            try:
                Sky_Arm_2 = np.array(fitsfiles[i+1][k].data)
                Sci_Arm_2 = np.array(fitsfiles[i][k].data)
                Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
                Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
                Sci_Arm_sub1 = []
                Sci_Arm_sub2 = []
                Sci_Arm_sub3 = []
                Sci_Arm_sub4 = []
                for j in range(len(Sky_Arm_2)):
                    Sci_Arm_sub1.append([])
                    Sci_Arm_sub2.append([])
                    Sci_Arm_sub3.append([])
                    Sci_Arm_sub4.append([])
                    new_flux1 = Sci_Arm_2[j] - f1_list[j] * Sky_Arm_2[j]
                    new_flux2 = Sci_Arm_2[j] - f2_list[j] * Sky_Arm_2[j]
                    new_flux3 = Sci_Arm_2[j] - f3_list[j] * Sky_Arm_2[j]
                    new_flux4 = Sci_Arm_2[j] - f4_list[j] * Sky_Arm_2[j]
                    Sci_Arm_sub1[j] = new_flux1
                    Sci_Arm_sub2[j] = new_flux2
                    Sci_Arm_sub3[j] = new_flux3
                    Sci_Arm_sub4[j] = new_flux4
                sub1 = np.mean(Sci_Arm_sub1)
                sub2 = np.mean(Sci_Arm_sub2)
                sub3 = np.mean(Sci_Arm_sub3)
                sub4 = np.mean(Sci_Arm_sub4)
                if sub1 == min(sub1, sub2, sub3, sub4):
                    Sci_Arm_sub = Sci_Arm_sub1
                elif sub2 == min(sub1, sub2, sub3, sub4):
                    Sci_Arm_sub = Sci_Arm_sub2
                elif sub3 == min(sub1, sub2, sub3, sub4):
                    Sci_Arm_sub = Sci_Arm_sub3
                elif sub4 == min(sub1, sub2, sub3, sub4):
                    Sci_Arm_sub = Sci_Arm_sub4
                Sci_Arm_sub = np.array(Sci_Arm_sub)
                fitsfiles[i][k].data = Sci_Arm_sub
            except:
                print "no data for arm "+str(k)
            #hdu = fits.PrimaryHDU(header=fitsfiles[i][0].header)
            #hdulist_file = fits.HDUList([hdu])
            #hdulist_file.append(fitsfiles[i][k])
            #hdulist_file.writeto('/home/rburnet/reflex/project_data/sky_sub/'+hdulist[i][-66:-5]+'_arm'+str(k)+'.fits')
        fitsfiles[i].writeto('/home/rburnet/reflex/project_data/sky_sub/'+hdulist[i][-66:])

    #CL0034 BA pairs
    if i == 2 or i == 7 or i == 18 or i == 23:
        Sky_Arm_2 = np.array(fitsfiles[i-1][2].data)
        Sci_Arm_2 = np.array(fitsfiles[i][2].data)
        Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
        Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
        f1_list = []
        for k in range(len(Sky_Arm_2)):
            f1_list.append([])
            Sky_flux = np.mean(Sky_Arm_2[k])
            Sci_flux = np.mean(Sci_Arm_2[k])
            f = Sci_flux / Sky_flux
            f1_list[k]  = f

        if i < 15:
            Sky_Arm_2 = np.array(fitsfiles[i-1][5].data)
            Sci_Arm_2 = np.array(fitsfiles[i][5].data)
            Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
            Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
            f2_list = []
            for k in range(len(Sky_Arm_2)):
                f2_list.append([])
                Sky_flux = np.mean(Sky_Arm_2[k])
                Sci_flux = np.mean(Sci_Arm_2[k])
                f = Sci_flux / Sky_flux
                f2_list[k]  = f
        else:
            f2_list = f1_list

    
        for k in range(len(fitsfiles[i])):
            try:
                Sky_Arm_2 = np.array(fitsfiles[i-1][k].data)
                Sci_Arm_2 = np.array(fitsfiles[i][k].data)
                Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
                Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
                Sci_Arm_sub1 = []
                Sci_Arm_sub2 = []
                for j in range(len(Sky_Arm_2)):
                    Sci_Arm_sub1.append([])
                    Sci_Arm_sub2.append([])
                    new_flux1 = Sci_Arm_2[j] - f1_list[j] * Sky_Arm_2[j]
                    new_flux2 = Sci_Arm_2[j] - f2_list[j] * Sky_Arm_2[j]
                    Sci_Arm_sub1[j] = new_flux1
                    Sci_Arm_sub2[j] = new_flux2
                sub1 = np.mean(Sci_Arm_sub1)
                sub2 = np.mean(Sci_Arm_sub2)
                if sub1 > sub2:
                    Sci_Arm_sub = Sci_Arm_sub1
                else:
                    Sci_Arm_sub = Sci_Arm_sub2
                Sci_Arm_sub = np.array(Sci_Arm_sub)
                fitsfiles[i][k].data = Sci_Arm_sub
            except:
                print "no data for arm "+str(k)
            #hdu = fits.PrimaryHDU(header=fitsfiles[i][0].header)
            #hdulist_file = fits.HDUList([hdu])
            #hdulist_file.append(fitsfiles[i][k])
            #hdulist_file.writeto('/home/rburnet/reflex/project_data/sky_sub/'+hdulist[i][-66:-5]+'_arm'+str(k)+'.fits')
        fitsfiles[i].writeto('/home/rburnet/reflex/project_data/sky_sub/'+hdulist[i][-66:])

    #CL0036 BA pairs
    if i == 10 or i == 15 or i == 28 or i == 33:
        if i < 16:
            Sky_Arm_2 = np.array(fitsfiles[i-1][5].data)
            Sci_Arm_2 = np.array(fitsfiles[i][5].data)
        if i > 16:
            Sky_Arm_2 = np.array(fitsfiles[i-1][9].data)
            Sci_Arm_2 = np.array(fitsfiles[i][9].data)
        Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
        Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
        f1_list = []
        for k in range(len(Sky_Arm_2)):
            f1_list.append([])
            Sky_flux = np.mean(Sky_Arm_2[k])
            Sci_flux = np.mean(Sci_Arm_2[k])
            f = Sci_flux / Sky_flux
            f1_list[k]  = f
        
        if i < 16:
            Sky_Arm_2 = np.array(fitsfiles[i-1][17].data)
            Sci_Arm_2 = np.array(fitsfiles[i][17].data)
        if i > 16:
            Sky_Arm_2 = np.array(fitsfiles[i-1][13].data)
            Sci_Arm_2 = np.array(fitsfiles[i][13].data)
        Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
        Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
        f2_list = []
        for k in range(len(Sky_Arm_2)):
            f2_list.append([])
            Sky_flux = np.mean(Sky_Arm_2[k])
            Sci_flux = np.mean(Sci_Arm_2[k])
            f = Sci_flux / Sky_flux
            f2_list[k]  = f

        Sky_Arm_2 = np.array(fitsfiles[i-1][21].data)
        Sci_Arm_2 = np.array(fitsfiles[i][21].data)
        Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
        Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
        f3_list = []
        for k in range(len(Sky_Arm_2)):
            f3_list.append([])
            Sky_flux = np.mean(Sky_Arm_2[k])
            Sci_flux = np.mean(Sci_Arm_2[k])
            f = Sci_flux / Sky_flux
            f3_list[k]  = f

        if i < 16:
            Sky_Arm_2 = np.array(fitsfiles[i-1][13].data)
            Sci_Arm_2 = np.array(fitsfiles[i][13].data)
            Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
            Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
            f4_list = []
            for k in range(len(Sky_Arm_2)):
                f4_list.append([])
                Sky_flux = np.mean(Sky_Arm_2[k])
                Sci_flux = np.mean(Sci_Arm_2[k])
                f = Sci_flux / Sky_flux
                f4_list[k]  = f
        else:
            f4_list = f1_list


        for k in range(len(fitsfiles[i])):
            try:
                Sky_Arm_2 = np.array(fitsfiles[i-1][k].data)
                Sci_Arm_2 = np.array(fitsfiles[i][k].data)
                Sky_Arm_2 = np.nan_to_num(Sky_Arm_2)
                Sci_Arm_2 = np.nan_to_num(Sci_Arm_2)
                Sci_Arm_sub1 = []
                Sci_Arm_sub2 = []
                Sci_Arm_sub3 = []
                Sci_Arm_sub4 = []
                for j in range(len(Sky_Arm_2)):
                    Sci_Arm_sub1.append([])
                    Sci_Arm_sub2.append([])
                    Sci_Arm_sub3.append([])
                    Sci_Arm_sub4.append([])
                    new_flux1 = Sci_Arm_2[j] - f1_list[j] * Sky_Arm_2[j]
                    new_flux2 = Sci_Arm_2[j] - f2_list[j] * Sky_Arm_2[j]
                    new_flux3 = Sci_Arm_2[j] - f3_list[j] * Sky_Arm_2[j]
                    new_flux4 = Sci_Arm_2[j] - f4_list[j] * Sky_Arm_2[j]
                    Sci_Arm_sub1[j] = new_flux1
                    Sci_Arm_sub2[j] = new_flux2
                    Sci_Arm_sub3[j] = new_flux3
                    Sci_Arm_sub4[j] = new_flux4
                sub1 = np.mean(Sci_Arm_sub1)
                sub2 = np.mean(Sci_Arm_sub2)
                sub3 = np.mean(Sci_Arm_sub3)
                sub4 = np.mean(Sci_Arm_sub4)
                if sub1 == min(sub1, sub2, sub3, sub4):
                    Sci_Arm_sub = Sci_Arm_sub1
                elif sub2 == min(sub1, sub2, sub3, sub4):
                    Sci_Arm_sub = Sci_Arm_sub2
                elif sub3 == min(sub1, sub2, sub3, sub4):
                    Sci_Arm_sub = Sci_Arm_sub3
                elif sub4 == min(sub1, sub2, sub3, sub4):
                    Sci_Arm_sub = Sci_Arm_sub4
                Sci_Arm_sub = np.array(Sci_Arm_sub)
                fitsfiles[i][k].data = Sci_Arm_sub
            except:
                    print "no data for arm "+str(k)
            #hdu = fits.PrimaryHDU(header=fitsfiles[i][0].header)
            #hdulist_file = fits.HDUList([hdu])
            #hdulist_file.append(fitsfiles[i][k])
            #hdulist_file.writeto('/home/rburnet/reflex/project_data/sky_sub/'+hdulist[i][-66:-5]+'_arm'+str(k)+'.fits')
        fitsfiles[i].writeto('/home/rburnet/reflex/project_data/sky_sub/'+hdulist[i][-66:])

