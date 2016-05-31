#!/usr/bin/python

####
# script to calculate f = <T_lambda(target)>/<S_lambda(sky)> , the fraction of average flux of the target frame target exposures to the average flux of the associated sky frame sky exposures to see if they agree. Then, it carries out sky subtraction of target cubes manually using T_lambda - S_lambda*f for each spaxel in the target cube and plots the spectrum of that and the scaled average sky spectrum in the range of +/- 0.05 microns from Halpha and saves them to ../figures/laptop/target_sky_subtraction_of_object/sky_subtracted_spectrum_of_bright_object_0.05microns/.
####

import numpy as np
import scipy as sp
import glob
import os
import fnmatch
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy.io import fits

hdulist_sci = []    #target cubes list
hdulist_sky = []    #associated sky cubes list

#append non-sky subtraction and non-sky tweaked target cubes to hdulist_sci list
for root,dirnames,filenames in os.walk('/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/'):
    for filename in fnmatch.filter(filenames, 'CL003[0-9]-YJ-OB-[1|2]_COMBINE_SCI_RECONSTRUCTED_[0-9]*.fits'):
        hdulist_sci.append(os.path.join(root, filename))

#append associated sky cubes for each target cube above to hdulist_sky list
for root,dirnames,filenames in os.walk('/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/'):
    for filename in fnmatch.filter(filenames, 'CL003[0-9]-YJ-OB-[1|2]_COMBINE_SCI_RECONSTRUCTED_ARM*_SKY.fits'):
        if not os.path.isfile(os.path.join(root, filename[:-8]+'SCI.fits')): #must check that sky cube doesn't have an associated sci cube, otherwise it is a sky cube for a sky arm in the science frame, not a sky cube for a target arm in the science frame.
            hdulist_sky.append(os.path.join(root, filename))

#really ugly method to associate sky arm names to target names by making a list of tuples.
for i in range(len(hdulist_sky)):
    if 'CL0034-YJ-OB-1' in hdulist_sky[i]:
        if 'ARM1_' in hdulist_sky[i]:
            hdulist_sky[i] = (1102, hdulist_sky[i])
        if 'ARM2_' in hdulist_sky[i]:
            hdulist_sky[i] = (1110, hdulist_sky[i])
        if 'ARM3_' in hdulist_sky[i]:
            hdulist_sky[i] = (1104, hdulist_sky[i])
        if 'ARM5_' in hdulist_sky[i]:
            hdulist_sky[i] = (1108, hdulist_sky[i])
        if 'ARM8_' in hdulist_sky[i]:
            hdulist_sky[i] = (1114, hdulist_sky[i])
        if 'ARM9_' in hdulist_sky[i]:
            hdulist_sky[i] = (1126, hdulist_sky[i])
        if 'ARM10_' in hdulist_sky[i]:
            hdulist_sky[i] = (1105, hdulist_sky[i])
        if 'ARM11_' in hdulist_sky[i]:
            hdulist_sky[i] = (1122, hdulist_sky[i])
        if 'ARM12_' in hdulist_sky[i]:
            hdulist_sky[i] = (1121, hdulist_sky[i])
        if 'ARM14_' in hdulist_sky[i]:
            hdulist_sky[i] = (1120, hdulist_sky[i])
        if 'ARM15_' in hdulist_sky[i]:
            hdulist_sky[i] = (1112, hdulist_sky[i])
        if 'ARM16_' in hdulist_sky[i]:
            hdulist_sky[i] = (1107, hdulist_sky[i])
        if 'ARM18_' in hdulist_sky[i]:
            hdulist_sky[i] = (1131, hdulist_sky[i])
        if 'ARM19_' in hdulist_sky[i]:
            hdulist_sky[i] = (1116, hdulist_sky[i])
        if 'ARM20_' in hdulist_sky[i]:
            hdulist_sky[i] = (1123, hdulist_sky[i])
        if 'ARM21_' in hdulist_sky[i]:
            hdulist_sky[i] = (1124, hdulist_sky[i])
        if 'ARM22_' in hdulist_sky[i]:
            hdulist_sky[i] = (1111, hdulist_sky[i])
        if 'ARM23_' in hdulist_sky[i]:
            hdulist_sky[i] = (1125, hdulist_sky[i])
        if 'ARM24_' in hdulist_sky[i]:
            hdulist_sky[i] = (1106, hdulist_sky[i])
    if 'CL0036-YJ-OB-1' in hdulist_sky[i]:
        if 'ARM1_' in hdulist_sky[i]:
            hdulist_sky[i] = (1230, hdulist_sky[i])
        if 'ARM2_' in hdulist_sky[i]:
            hdulist_sky[i] = (1204, hdulist_sky[i])
        if 'ARM3_' in hdulist_sky[i]:
            hdulist_sky[i] = (1217, hdulist_sky[i])
        if 'ARM4_' in hdulist_sky[i]:
            hdulist_sky[i] = (1212, hdulist_sky[i])
        if 'ARM5_' in hdulist_sky[i]:
            hdulist_sky[i] = (1228, hdulist_sky[i])
        if 'ARM8_' in hdulist_sky[i]:
            hdulist_sky[i] = (1206, hdulist_sky[i])
        if 'ARM9_' in hdulist_sky[i]:
            hdulist_sky[i] = (1219, hdulist_sky[i])
        if 'ARM10_' in hdulist_sky[i]:
            hdulist_sky[i] = (1211, hdulist_sky[i])
        if 'ARM11_' in hdulist_sky[i]:
            hdulist_sky[i] = (1215, hdulist_sky[i])
        if 'ARM12_' in hdulist_sky[i]:
            hdulist_sky[i] = (1203, hdulist_sky[i])
        if 'ARM14_' in hdulist_sky[i]:
            hdulist_sky[i] = (1234, hdulist_sky[i])
        if 'ARM16_' in hdulist_sky[i]:
            hdulist_sky[i] = (1208, hdulist_sky[i])
        if 'ARM18_' in hdulist_sky[i]:
            hdulist_sky[i] = (1205, hdulist_sky[i])
        if 'ARM20_' in hdulist_sky[i]:
            hdulist_sky[i] = (1210, hdulist_sky[i])
        if 'ARM21_' in hdulist_sky[i]:
            hdulist_sky[i] = (1218, hdulist_sky[i])
        if 'ARM22_' in hdulist_sky[i]:
            hdulist_sky[i] = (1207, hdulist_sky[i])
        if 'ARM24_' in hdulist_sky[i]:
            hdulist_sky[i] = (1213, hdulist_sky[i])
    if 'CL0034-YJ-OB-2' in hdulist_sky[i]:
        if 'ARM1_' in hdulist_sky[i]:
            hdulist_sky[i] = (2102, hdulist_sky[i])
        if 'ARM2_' in hdulist_sky[i]:
            hdulist_sky[i] = (2110, hdulist_sky[i])
        if 'ARM3_' in hdulist_sky[i]:
            hdulist_sky[i] = (2104, hdulist_sky[i])
        if 'ARM5_' in hdulist_sky[i]:
            hdulist_sky[i] = (2108, hdulist_sky[i])
        if 'ARM8_' in hdulist_sky[i]:
            hdulist_sky[i] = (2114, hdulist_sky[i])
        if 'ARM9_' in hdulist_sky[i]:
            hdulist_sky[i] = (2126, hdulist_sky[i])
        if 'ARM10_' in hdulist_sky[i]:
            hdulist_sky[i] = (2105, hdulist_sky[i])
        if 'ARM11_' in hdulist_sky[i]:
            hdulist_sky[i] = (2122, hdulist_sky[i])
        if 'ARM12_' in hdulist_sky[i]:
            hdulist_sky[i] = (2121, hdulist_sky[i])
        if 'ARM14_' in hdulist_sky[i]:
            hdulist_sky[i] = (2120, hdulist_sky[i])
        if 'ARM15_' in hdulist_sky[i]:
            hdulist_sky[i] = (2112, hdulist_sky[i])
        if 'ARM16_' in hdulist_sky[i]:
            hdulist_sky[i] = (2107, hdulist_sky[i])
        if 'ARM18_' in hdulist_sky[i]:
            hdulist_sky[i] = (2131, hdulist_sky[i])
        if 'ARM19_' in hdulist_sky[i]:
            hdulist_sky[i] = (2116, hdulist_sky[i])
        if 'ARM20_' in hdulist_sky[i]:
            hdulist_sky[i] = (2123, hdulist_sky[i])
        if 'ARM21_' in hdulist_sky[i]:
            hdulist_sky[i] = (2124, hdulist_sky[i])
        if 'ARM22_' in hdulist_sky[i]:
            hdulist_sky[i] = (2111, hdulist_sky[i])
        if 'ARM23_' in hdulist_sky[i]:
            hdulist_sky[i] = (2125, hdulist_sky[i])
        if 'ARM24_' in hdulist_sky[i]:
            hdulist_sky[i] = (2106, hdulist_sky[i])
    if 'CL0036-YJ-OB-2' in hdulist_sky[i]:
        if 'ARM1_' in hdulist_sky[i]:
            hdulist_sky[i] = (2230, hdulist_sky[i])
        if 'ARM2_' in hdulist_sky[i]:
            hdulist_sky[i] = (2204, hdulist_sky[i])
        if 'ARM3_' in hdulist_sky[i]:
            hdulist_sky[i] = (2217, hdulist_sky[i])
        if 'ARM4_' in hdulist_sky[i]:
            hdulist_sky[i] = (2212, hdulist_sky[i])
        if 'ARM5_' in hdulist_sky[i]:
            hdulist_sky[i] = (2228, hdulist_sky[i])
        if 'ARM8_' in hdulist_sky[i]:
            hdulist_sky[i] = (2206, hdulist_sky[i])
        if 'ARM9_' in hdulist_sky[i]:
            hdulist_sky[i] = (2219, hdulist_sky[i])
        if 'ARM10_' in hdulist_sky[i]:
            hdulist_sky[i] = (2211, hdulist_sky[i])
        if 'ARM11_' in hdulist_sky[i]:
            hdulist_sky[i] = (2215, hdulist_sky[i])
        if 'ARM12_' in hdulist_sky[i]:
            hdulist_sky[i] = (2203, hdulist_sky[i])
        if 'ARM14_' in hdulist_sky[i]:
            hdulist_sky[i] = (2234, hdulist_sky[i])
        if 'ARM16_' in hdulist_sky[i]:
            hdulist_sky[i] = (2208, hdulist_sky[i])
        if 'ARM18_' in hdulist_sky[i]:
            hdulist_sky[i] = (2205, hdulist_sky[i])
        if 'ARM20_' in hdulist_sky[i]:
            hdulist_sky[i] = (2210, hdulist_sky[i])
        if 'ARM21_' in hdulist_sky[i]:
            hdulist_sky[i] = (2218, hdulist_sky[i])
        if 'ARM22_' in hdulist_sky[i]:
            hdulist_sky[i] = (2207, hdulist_sky[i])
        if 'ARM24_' in hdulist_sky[i]:
            hdulist_sky[i] = (2213, hdulist_sky[i])

#now I need to sort hdulist_sci with the same indexes. This is thankfully a little cleaner.
for i in range(len(hdulist_sci)):
    if 'CL0034-YJ-OB-1' in hdulist_sci[i]:
        if '_' in hdulist_sci[i][-7:-5]:
            hdulist_sci[i] = (1100+int(hdulist_sci[i][-6:-5]), hdulist_sci[i])
        else:
            hdulist_sci[i] = (1100+int(hdulist_sci[i][-7:-5]), hdulist_sci[i])
    if 'CL0036-YJ-OB-1' in hdulist_sci[i]:
        if '_' in hdulist_sci[i][-7:-5]:
            hdulist_sci[i] = (1200+int(hdulist_sci[i][-6:-5]), hdulist_sci[i])
        else:
            hdulist_sci[i] = (1200+int(hdulist_sci[i][-7:-5]), hdulist_sci[i]) 
    if 'CL0034-YJ-OB-2' in hdulist_sci[i]:
        if '_' in hdulist_sci[i][-7:-5]:
            hdulist_sci[i] = (2100+int(hdulist_sci[i][-6:-5]), hdulist_sci[i])
        else:
            hdulist_sci[i] = (2100+int(hdulist_sci[i][-7:-5]), hdulist_sci[i])
    if 'CL0036-YJ-OB-2' in hdulist_sci[i]:
        if '_' in hdulist_sci[i][-7:-5]:
            hdulist_sci[i] = (2200+int(hdulist_sci[i][-6:-5]), hdulist_sci[i])
        else:
            hdulist_sci[i] = (2200+int(hdulist_sci[i][-7:-5]), hdulist_sci[i])

#Where the sorting happens
hdulist_sci = sorted(hdulist_sci)
hdulist_sky = sorted(hdulist_sky)

#this will return the lists as lists of filenames, but now they are in the right order. ie. hdulist_sci[1] is the target cube and it's associated sky cube is hdulist_sky[1].
for i in range(len(hdulist_sci)):
    hdulist_sci[i] = hdulist_sci[i][1]
    hdulist_sky[i] = hdulist_sky[i][1]

#Arrays to select perimeter spaxels only from target cubes, to be used in the fraction.
spaxel_array_select_OB1 = []
spaxel_array_select_OB2 = []
spaxel_array_select_OB2_special = []

for i in range(14):
    if i == 0 or i == 13:
        spaxel_array_select_OB1.append([False, False, False, False, False, False, False, False, False, False, False, False, False, False, False, False])
    elif i < 3 or i > 10:
        spaxel_array_select_OB1.append([False, True, True, True, True, True, True, True, True, True, True, True, True, True, True, False])
    else:
        spaxel_array_select_OB1.append([False, True, True, False, False, False, False, False, False, False, False, False, False, True, True, False])

for i in range(15):
    if i == 0 or i == 14:
        spaxel_array_select_OB2.append([False, False, False, False, False, False, False, False, False, False, False, False, False, False])
    elif i < 3 or i > 11:
        spaxel_array_select_OB2.append([False, True, True, True, True, True, True, True, True, True, True, True, True, False])
    else:
        spaxel_array_select_OB2.append([False, True, True, False, False, False, False, False, False, False, False, True, True, False])

for i in range(14):
    if i == 0 or i == 13:
        spaxel_array_select_OB2_special.append([False, False, False, False, False, False, False, False, False, False, False, False, False, False])
    elif i < 3 or i > 10:
        spaxel_array_select_OB2_special.append([False, True, True, True, True, True, True, True, True, True, True, True, True, False])
    else:
        spaxel_array_select_OB2_special.append([False, True, True, False, False, False, False, False, False, False, False, True, True, False])

spaxel_array_select_OB1 = np.array(spaxel_array_select_OB1)
spaxel_array_select_OB2 = np.array(spaxel_array_select_OB2)
spaxel_array_select_OB2_special = np.array(spaxel_array_select_OB2_special)

f = []  #inital fraction for each wavelength slice of every pair of files. Has shape (#of pairs of files, 2048)

for i in range(len(hdulist_sci)):
    f.append([])
    hdufits_sci = fits.open(hdulist_sci[i])
    hdufits_sky = fits.open(hdulist_sky[i])
    print hdulist_sci[i]
    for j in range(len(hdufits_sci[1].data)):
        spaxel_array = np.nan_to_num(hdufits_sci[1].data[j])
        if 'OB-1' in hdulist_sci[i]:
            spaxel_array = spaxel_array[spaxel_array_select_OB1]
        elif 'CL0036-YJ-OB-2' in hdulist_sci[i]:
            spaxel_array = spaxel_array[spaxel_array_select_OB2_special]
        else:
            spaxel_array = spaxel_array[spaxel_array_select_OB2]
        avg_flux_sci = np.mean(spaxel_array)
        #avg_flux_sci = np.mean(hdufits_sci[1].data[j][~np.isnan(hdufits_sci[1].data[j])])
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

#Now, let's do sky subtraction of target cube manually by carrying out T_lambda - S_lambda*f for every spaxel in the T_lambda (target) cube and then plot the spectrum.
sky_subtracted_sum_data = []    #sum of sky subtracted data for each wavelength slice, used to plot spectrum
avg_flux_sky_list = []          #average sky flux list, used to plot sky flux spectrum

for i in range(len(hdulist_sci)):
    sky_subtracted_sum_data.append([])
    avg_flux_sky_list.append([])
    hdufits_sci = fits.open(hdulist_sci[i])
    hdufits_sky = fits.open(hdulist_sky[i])
    for j in range(len(hdufits_sci[1].data)):
        avg_flux_sky = np.mean(hdufits_sky[1].data[j][~np.isnan(hdufits_sky[1].data[j])])
        spaxel_array = np.array(hdufits_sci[1].data[j])
        spaxel_array = spaxel_array - avg_flux_sky*f[i][j]
   #     total_flux_sky = np.sum(hdufits_sky[1].data[j][~np.isnan(hdufits_sky[1].data[j])])
        sky_subtracted_sum_data[i].append(np.sum(spaxel_array[~np.isnan(spaxel_array)]))
        avg_flux_sky_list[i].append(avg_flux_sky*f[i][j])
   #     hdufits_sci[1].data[j] = spaxel_array
   # hdufits_sci.writeto('/home/rburnet/reflex/project_data/test_cubes/'+hdulist_sci[i][-48:])   #write sky subtracted cubes
    hdufits_sci.close()
    hdufits_sky.close()

for i in range(len(sky_subtracted_sum_data)):
    if 'CL0034-YJ-OB-1_COMBINE_SCI_RECONSTRUCTED_11' in hdulist_sci[i] or 'CL0034-YJ-OB-1_COMBINE_SCI_RECONSTRUCTED_5' in hdulist_sci[i] or 'CL0036-YJ-OB-1_COMBINE_SCI_RECONSTRUCTED_4' in hdulist_sci[i] or 'CL0036-YJ-OB-1_COMBINE_SCI_RECONSTRUCTED_15' in hdulist_sci[i]:
        #loop to turn all NaN values in sky_subtracted_sum_data array to 0.0
        for j in range(len(sky_subtracted_sum_data[i])):
            if sky_subtracted_sum_data[i][j] != sky_subtracted_sum_data[i][j]:
                sky_subtracted_sum_data[i][j] = 0.0
        if not os.path.exists('../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_bright_object_0.05microns/'):
            os.makedirs('../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_bright_object_0.05microns/')
        line1, = plt.plot(x,sky_subtracted_sum_data[i], label='Total Flux of Target Slice')
        if 'RECONSTRUCTED_5' in hdulist_sci[i]:
            plt.title('Total Flux vs Wavelength  from \n '+hdulist_sci[i][-47:], fontsize=10)
            z = 0.869
        if 'RECONSTRUCTED_11' in hdulist_sci[i]:
            plt.title('Total Flux vs Wavelength  from \n '+hdulist_sci[i][-48:], fontsize=10)
            z = 0.867
        if 'RECONSTRUCTED_4' in hdulist_sci[i]:
            plt.title('Total Flux vs Wavelength  from \n '+hdulist_sci[i][-47:], fontsize=10)
            z = 0.861
        if 'RECONSTRUCTED_15' in hdulist_sci[i]:
            plt.title('Total Flux vs Wavelength  from \n '+hdulist_sci[i][-48:], fontsize=10)
            z = 0.868
        l = 0.6563 * ( 1 + z )
        NII_1 = 0.6548 * ( 1 + z )
        NII_2 = 0.6584 * ( 1 + z )
        line2, = plt.plot((l,l),(-3e-16,3e16))
        line3, = plt.plot((NII_1,NII_1),(-3e-16,3e16))
        line4, = plt.plot((NII_2,NII_2),(-3e-16,3e16))
        line5, = plt.plot((1.0,1.35882),(0.0,0.0))
        line6, = plt.plot(x,avg_flux_sky_list[i], label='Scaled Average Sky Flux')
        plt.xlabel('Wavelength (microns)',fontsize=10)
        plt.ylabel('Total Flux (ergs cm$\mathregular{^{-2}}$ $\\AA\mathregular{^{-1}}$ s$\mathregular{^{-1}}$)',fontsize=10)
        plt.axis('tight')
        plt.ylim([-3e-16,3e-16])
        plt.xlim([l-0.05,l+0.05])
        fig = plt.gcf()
        lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2)
        fig.set_size_inches(18.5, 5.5)
        plt.setp(line1, linewidth=0.1, color='b')
        plt.setp(line2, linewidth=0.1)
        plt.setp(line3, linewidth=0.1)
        plt.setp(line4, linewidth=0.1)
        plt.setp(line5, linewidth=0.1)
        plt.setp(line6, linewidth=0.1)
        print '../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_bright_object_0.05microns/'+hdulist_sci[i][-48:]+'.pdf'
        plt.savefig('../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_bright_object_0.05microns/'+hdulist_sci[i][-48:]+'.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')
        plt.close()
