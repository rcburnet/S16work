#!/usr/bin/python

####
# script to calculate f = <T_lambda(target)>/<S_lambda(sky)> , the fraction of average flux of the target frame target exposures to the average flux of the associated sky frame sky exposures to see if they agree on an IFU by IFU basis so that the cubes are the same shape ONLY FOR OB1. Then, it carries out sky subtraction of target cubes manually using T_lambda - S_lambda*f for each spaxel in the target cube and plots the mean spectrum of that and the scaled mean sky spectrum in the range of +/- 0.05 microns from Halpha and saves them to ../figures/laptop/target_sky_subtraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets/.
####

import numpy as np
import scipy as sp
import glob
import os
import fnmatch
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy.io import fits

hdulist_sci = [['/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-07-11T06:34:52.380_tpl/CL0034-YJ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-07-11T06:34:52.380.fits', '/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-07-11T06:34:52.380_tpl/CL0034-YJ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-07-11T06:50:52.037.fits', '/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-07-11T06:34:52.380_tpl/CL0034-YJ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-07-11T06:59:01.745.fits'],['/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-07-11T08:24:20.562_tpl/CL0036-YJ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-07-11T08:24:20.562.fits', '/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-07-11T08:24:20.562_tpl/CL0036-YJ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-07-11T08:40:15.105.fits', '/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-07-11T08:24:20.562_tpl/CL0036-YJ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-07-11T08:48:01.174.fits']]  #science list

hdulist_sky = [['/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-07-11T06:34:52.380_tpl/CL0034-YJ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-07-11T06:42:39.007.fits', '/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-07-11T06:34:52.380_tpl/CL0034-YJ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-07-11T07:06:49.051.fits'],['/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-07-11T08:24:20.562_tpl/CL0036-YJ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-07-11T08:32:06.663.fits', '/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-07-11T08:24:20.562_tpl/CL0036-YJ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-07-11T08:55:48.356.fits']]   #associated sky frames list

spaxel_array_select_OB1 = []    #array to select perimeter pixels of target cubes

for i in range(14):
    if i == 0 or i == 13:
        spaxel_array_select_OB1.append([False, False, False, False, False, False, False, False, False, False, False,  False, False, False])
    elif i < 3 or i > 10:
        spaxel_array_select_OB1.append([False, True, True, True, True, True, True, True, True, True, True, True, True, False])
    else:
        spaxel_array_select_OB1.append([False, True, True, False, False, False, False, False, False, False, False, True, True, False])

spaxel_array_select_OB1 = np.array(spaxel_array_select_OB1)

x = np.linspace(1,1.35882,2048)

for i in range(len(hdulist_sci)):
    for j in range(len(hdulist_sci[i])):
        hdufits_sci = fits.open(hdulist_sci[i][j])  #open science frame
        if j == 0 or j == 1:    #first two science frames
            hdufits_sky = fits.open(hdulist_sky[i][0])  #open first sky frame
        if j == 2:              #last science frame
            hdufits_sky = fits.open(hdulist_sky[i][1])  #open last sky frame
        #for k in range(len(hdufits_sci)):   #loop into each ifu of the science/sky frame
        #    sky_subtracted_sum_data = []
        #    avg_flux_sky_list = []
        if i == 0:
            k_list = [8,20]  #IFU 8 = Target 5, IFU 20 = Target 11
            target = {8:5, 20:11}
        if i == 1:
            k_list = [24,9]  #IFU 24 = Target 4, IFU 9 = Target 15
            target = {24:4, 9:15}
        for k in k_list:
            sky_subtracted_sum_data = []
            avg_flux_sky_list = []
            if not type(hdufits_sci[k].data) == type(hdufits_sci[0].data):  #check to make sure there is a data array for that specific IFU
                for m in range(len(hdufits_sci[k].data)):
                    spaxel_array = np.nan_to_num(hdufits_sci[k].data[m])
                    spaxel_array = spaxel_array[spaxel_array_select_OB1]
                    avg_flux_sci = np.mean(spaxel_array)
                    avg_flux_sky = np.mean(hdufits_sky[k].data[m][~np.isnan(hdufits_sky[k].data[m])])
                    if abs(avg_flux_sci/avg_flux_sky) == np.inf:
                        f = np.nan
                    else:
                        f = avg_flux_sci/avg_flux_sky
                    if f != f:
                        f = 1.0
                    spaxel_array = np.array(hdufits_sci[k].data[m])
                    sky_array = np.array(hdufits_sky[k].data[m])
                    spaxel_array = spaxel_array - sky_array*f
                    sky_subtracted_sum_data.append(np.sum(spaxel_array[~np.isnan(spaxel_array)]))
                    avg_flux_sky_list.append(avg_flux_sky*f)
                    if sky_subtracted_sum_data[m] != sky_subtracted_sum_data[m]:
                        sky_subtracted_sum_data[m] = 0.0
                    if avg_flux_sky_list[m] != avg_flux_sky_list[m]:
                        avg_flux_sky_list[m] = 0.0
            #plot spectrum for each arm now
                if not os.path.exists('../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets/'):
                    os.makedirs('../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets/')
                if k == 8:
                    z = 0.869
                if k == 20:
                    z = 0.867
                if k == 24:
                    z = 0.861
                if k == 9:
                    z = 0.868
                y_range = (min(sky_subtracted_sum_data),max(sky_subtracted_sum_data))
                l = 0.6563 * ( 1 + z )
                NII_1 = 0.6548 * ( 1 + z )
                NII_2 = 0.6584 * ( 1 + z )
                line1, = plt.plot(x,sky_subtracted_sum_data, label='Total Flux of Target IFU')
                line2, = plt.plot((l,l),y_range)
                line3, = plt.plot((NII_1,NII_1),y_range)
                line4, = plt.plot((NII_2,NII_2),y_range)
                line5, = plt.plot((1.0,1.35882),(0.0,0.0))
                line6, = plt.plot(x, avg_flux_sky_list, label='Mean Flux of Scaled Sky IFU')
                plt.xlabel('Wavelength (microns)',fontsize=10)
                plt.ylabel('Flux (ergs cm$\mathregular{^{-2}}$ $\\AA\mathregular{^{-1}}$ s$\mathregular{^{-1}}$)',fontsize=10)
                plt.axis('tight')
                plt.ylim([-1e-16,1e-16])
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
                print '../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets/'+hdulist_sci[i][j][-66:-52]+'_'+hdulist_sci[i][j][-28:-5]+'_Target_'+str(target[k])+'.pdf'
                plt.savefig('../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets/'+hdulist_sci[i][j][-66:-52]+'_'+hdulist_sci[i][j][-28:-5]+'_Target_'+str(target[k])+'.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')
                plt.close()
        hdufits_sci.close()
        hdufits_sky.close()
        
