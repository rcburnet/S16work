#!/usr/bin/python

####
# script to calculate f = <T_lambda(target)>/<S_lambda(sky)> , the fraction of average flux of the target frame target exposures to the average flux of the associated sky frame sky exposures to see if they agree on an IFU by IFU basis so that the cubes are the same shape FOR BOTH OB1 AND OB2. Then, it carries out sky subtraction of target cubes manually using T_lambda - S_lambda*f for each spaxel in the target cube and plots the mean spectrum of that, averaged for all exposures, in the range of +/- 0.05 microns from Halpha and saves them to ../figures/laptop/target_sky_subtraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets_specific_regions_exposures_averaged/
####

import numpy as np
import scipy as sp
import glob
import os
import fnmatch
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy.io import fits

hdulist_sci = [['/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-08-12T03:46:35.963_tpl/CL0034-IZ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-08-12T03:46:35.963.fits', '/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-08-12T03:46:35.963_tpl/CL0034-IZ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-08-12T04:02:31.804.fits', '/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-08-12T03:46:35.963_tpl/CL0034-IZ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-08-12T04:10:14.875.fits'],['/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-08-13T03:45:05.815_tpl/CL0034-IZ-OB-2_SCI_RECONSTRUCTED_KMOS.2014-08-13T03:45:05.815.fits','/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-08-13T03:45:05.815_tpl/CL0034-IZ-OB-2_SCI_RECONSTRUCTED_KMOS.2014-08-13T04:01:01.929.fits','/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-08-13T03:45:05.815_tpl/CL0034-IZ-OB-2_SCI_RECONSTRUCTED_KMOS.2014-08-13T04:09:09.077.fits'],['/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-09-07T02:46:11.187_tpl/CL0036-IZ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-09-07T02:46:11.187.fits', '/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-09-07T02:46:11.187_tpl/CL0036-IZ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-09-07T03:02:02.009.fits', '/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-09-07T02:46:11.187_tpl/CL0036-IZ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-09-07T03:09:44.376.fits'],['/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-09-07T03:50:42.916_tpl/CL0036-IZ-OB-2_SCI_RECONSTRUCTED_KMOS.2014-09-07T03:50:42.916.fits','/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-09-07T03:50:42.916_tpl/CL0036-IZ-OB-2_SCI_RECONSTRUCTED_KMOS.2014-09-07T04:06:32.653.fits','/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-09-07T03:50:42.916_tpl/CL0036-IZ-OB-2_SCI_RECONSTRUCTED_KMOS.2014-09-07T04:14:15.736.fits']]  #science list

hdulist_sky = [['/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-08-12T03:46:35.963_tpl/CL0034-IZ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-08-12T03:54:21.205.fits', '/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-08-12T03:46:35.963_tpl/CL0034-IZ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-08-12T04:18:00.831.fits'],['/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-08-13T03:45:05.815_tpl/CL0034-IZ-OB-2_SCI_RECONSTRUCTED_KMOS.2014-08-13T03:52:51.739.fits','/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-08-13T03:45:05.815_tpl/CL0034-IZ-OB-2_SCI_RECONSTRUCTED_KMOS.2014-08-13T04:16:53.533.fits'],['/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-09-07T02:46:11.187_tpl/CL0036-IZ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-09-07T02:53:55.771.fits', '/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-09-07T02:46:11.187_tpl/CL0036-IZ-OB-1_SCI_RECONSTRUCTED_KMOS.2014-09-07T03:17:29.031.fits'],['/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-09-07T03:50:42.916_tpl/CL0036-IZ-OB-2_SCI_RECONSTRUCTED_KMOS.2014-09-07T03:58:27.536.fits','/home/rburnet/reflex/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/KMOS.2014-09-07T03:50:42.916_tpl/CL0036-IZ-OB-2_SCI_RECONSTRUCTED_KMOS.2014-09-07T04:22:01.121.fits']]   #associated sky frames list

spaxel_array_select_OB1 = []    #array to select perimeter pixels of target cubes
spaxel_array_select_target_11 = []
spaxel_array_select_target_5 = []
spaxel_array_select_target_4 = []

for i in range(14):
    if i == 0 or i == 13:
        spaxel_array_select_OB1.append([False, False, False, False, False, False, False, False, False, False, False,  False, False, False])
    elif i < 3 or i > 10:
        spaxel_array_select_OB1.append([False, True, True, True, True, True, True, True, True, True, True, True, True, False])
    else:
        spaxel_array_select_OB1.append([False, True, True, False, False, False, False, False, False, False, False, True, True, False])

spaxel_array_select_OB1 = np.array(spaxel_array_select_OB1)

for i in range(14):
    if i < 6 or i > 9:
        spaxel_array_select_target_11.append([False, False, False, False, False, False, False, False, False, False, False, False, False, False])
    else:
        spaxel_array_select_target_11.append([False, False, False, False, False, False, False, True, True, True, True, False, False, False])

spaxel_array_select_target_11 = np.array(spaxel_array_select_target_11)

for i in range(14):
    if i < 3 or i > 10:
        spaxel_array_select_target_5.append([False, False, False, False, False, False, False, False, False, False, False, False, False, False])
    else:
        spaxel_array_select_target_5.append([False, False, False, False, False, True, True, True, True, False, False, False, False, False])

spaxel_array_select_target_5 = np.array(spaxel_array_select_target_5)

for i in range(14):
    if i < 6 or i > 10:
        spaxel_array_select_target_4.append([False, False, False, False, False, False, False, False, False, False, False, False, False, False])
    else:
        spaxel_array_select_target_4.append([False, False, False, False, False, False, False, False, True, True, True, True, False, False])

spaxel_array_select_target_4 = np.array(spaxel_array_select_target_4)



x = np.linspace(0.78,1.08985,2048)

for i in range(len(hdulist_sci)):
    if i == 0 or i == 2:
        sky_subtracted_sum_data_tot1 = []
        sky_subtracted_sum_data_tot2 = []
    for j in range(len(hdulist_sci[i])):
        hdufits_sci = fits.open(hdulist_sci[i][j])  #open science frame
        if j == 0 or j == 1:    #first two science frames
            hdufits_sky = fits.open(hdulist_sky[i][0])  #open first sky frame
        if j == 2:              #last science frame
            hdufits_sky = fits.open(hdulist_sky[i][1])  #open last sky frame
        #if i == 1 or i == 3:
        #    hdufits_sky = fits.open(hdulist_sky[i][0])
        #for k in range(len(hdufits_sci)):   #loop into each ifu of the science/sky frame
        #    sky_subtracted_sum_data = []
        #    avg_flux_sky_list = []
        if i == 0 or i == 1:
            k_list = [8,20]  #IFU 8 = Target 5, IFU 20 = Target 11
            target = {8:5, 20:11}
        if i == 2 or i == 3:
            k_list = [16,1]  #IFU 24 = Target 4, IFU 9 = Target 15
            target = {16:4, 1:15}
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
                    spaxel_array = np.nan_to_num(spaxel_array)
                    if k == 8:
                        spaxel_array = spaxel_array#[spaxel_array_select_target_5]
                    if k == 20:
                        spaxel_array = spaxel_array#[spaxel_array_select_target_11]
                    if k == 16:
                        spaxel_array = spaxel_array#[spaxel_array_select_target_4]
                    if k == 1:
                        spaxel_array = spaxel_array#[spaxel_array_select_target_11]
                    sky_subtracted_sum_data.append(np.sum(spaxel_array[~np.isnan(spaxel_array)]))
                    if k == 8 or k == 16:
                        sky_subtracted_sum_data_tot1.append(sky_subtracted_sum_data)
                    else:
                        sky_subtracted_sum_data_tot2.append(sky_subtracted_sum_data)
                    avg_flux_sky_list.append(avg_flux_sky*f)
                    if sky_subtracted_sum_data[m] != sky_subtracted_sum_data[m]:
                        sky_subtracted_sum_data[m] = 0.0
                    if avg_flux_sky_list[m] != avg_flux_sky_list[m]:
                        avg_flux_sky_list[m] = 0.0
            #plot spectrum for each arm now
        hdufits_sci.close()
        hdufits_sky.close()
    if i == 1 or i == 3:
        sky_subtracted_sum_data_tot1 = np.mean(sky_subtracted_sum_data_tot1, axis=0)
        sky_subtracted_sum_data_tot2 = np.mean(sky_subtracted_sum_data_tot2, axis=0)
    if i == 1:
        if not os.path.exists('../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets_specific_regions_exposures_averaged_IZ/'):
            os.makedirs('../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets_specific_regions_exposures_averaged_IZ/')
        z = 0.869
        #if k == 20:
        #    z = 0.867
        #if k == 24:
        #    z = 0.861
        #if k == 9:
        #    z = 0.868
        y_range = (min(sky_subtracted_sum_data_tot1),max(sky_subtracted_sum_data_tot1))
        Hbeta = 0.4861 * ( 1 + z )
        OIII = 0.5007 * ( 1 + z )
        line1, = plt.plot(x,sky_subtracted_sum_data_tot1, label='Total Flux of Target IFU')
        line2, = plt.plot((Hbeta,Hbeta),y_range)
        line3, = plt.plot((OIII,OIII),y_range)
        line5, = plt.plot((0.78,1.08985),(0.0,0.0))
        #line6, = plt.plot(x, avg_flux_sky_list, label='Mean Flux of Scaled Sky IFU')
        plt.xlabel('Wavelength (microns)',fontsize=10)
        plt.ylabel('Flux (ergs cm$\mathregular{^{-2}}$ $\\AA\mathregular{^{-1}}$ s$\mathregular{^{-1}}$)',fontsize=10)
        plt.axis('tight')
        plt.ylim([-5e-16,5e-16])
        #plt.xlim([(Hbeta+OIII)/2.0-0.05,(Hbeta+OIII)/2.0+0.05])
        fig = plt.gcf()
        lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2)
        fig.set_size_inches(18.5, 5.5)
        plt.setp(line1, linewidth=0.1, color='b')
        plt.setp(line2, linewidth=0.1)
        plt.setp(line3, linewidth=0.1)
        plt.setp(line5, linewidth=0.1)
        #plt.setp(line6, linewidth=0.1)
        print '../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets_specific_regions_exposures_averaged_IZ/'+hdulist_sci[i][j][-66:-57]+'_Target_5_specific_region.pdf'
        plt.savefig('../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets_specific_regions_exposures_averaged_IZ/'+hdulist_sci[i][j][-66:-57]+'_Target_5_specific_region.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')

        plt.close()
        #z = 0.869
        z = 0.867
        #if k == 24:
        #    z = 0.861
        #if k == 9:
        #    z = 0.868
        y_range = (min(sky_subtracted_sum_data_tot2),max(sky_subtracted_sum_data_tot2))
        Hbeta = 0.4861 * ( 1 + z )
        OIII = 0.5007 * ( 1 + z )
        line1, = plt.plot(x,sky_subtracted_sum_data_tot2, label='Total Flux of Target IFU')
        line2, = plt.plot((Hbeta,Hbeta),y_range)
        line3, = plt.plot((OIII,OIII),y_range)
        line5, = plt.plot((0.78,1.08985),(0.0,0.0))
        #line6, = plt.plot(x, avg_flux_sky_list, label='Mean Flux of Scaled Sky IFU')
        plt.xlabel('Wavelength (microns)',fontsize=10)
        plt.ylabel('Flux (ergs cm$\mathregular{^{-2}}$ $\\AA\mathregular{^{-1}}$ s$\mathregular{^{-1}}$)',fontsize=10)
        plt.axis('tight')
        plt.ylim([-5e-16,5e-16])
        #plt.xlim([(Hbeta+OIII)/2.0-0.05,(Hbeta+OIII)/2.0+0.05])
        fig = plt.gcf()
        lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2)
        fig.set_size_inches(18.5, 5.5)
        plt.setp(line1, linewidth=0.1, color='b')
        plt.setp(line2, linewidth=0.1)
        plt.setp(line3, linewidth=0.1)
        plt.setp(line5, linewidth=0.1)
        #plt.setp(line6, linewidth=0.1)
        print '../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets_specific_regions_exposures_averaged_IZ/'+hdulist_sci[i][j][-66:-57]+'_Target_11_specific_region.pdf'
        plt.savefig('../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets_specific_regions_exposures_averaged_IZ/'+hdulist_sci[i][j][-66:-57]+'_Target_11_specific_region.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')

        plt.close()
    if i == 3:
        #z = 0.869
        #if k == 20:
        #    z = 0.867
        z = 0.861
        #if k == 9:
        #    z = 0.868
        y_range = (min(sky_subtracted_sum_data_tot1),max(sky_subtracted_sum_data_tot1))
        Hbeta = 0.4861 * ( 1 + z )
        OIII = 0.5007 * ( 1 + z )
        line1, = plt.plot(x,sky_subtracted_sum_data_tot1, label='Total Flux of Target IFU')
        line2, = plt.plot((Hbeta,Hbeta),y_range)
        line3, = plt.plot((OIII,OIII),y_range)
        line5, = plt.plot((0.78,1.08985),(0.0,0.0))
        #line6, = plt.plot(x, avg_flux_sky_list, label='Mean Flux of Scaled Sky IFU')
        plt.xlabel('Wavelength (microns)',fontsize=10)
        plt.ylabel('Flux (ergs cm$\mathregular{^{-2}}$ $\\AA\mathregular{^{-1}}$ s$\mathregular{^{-1}}$)',fontsize=10)
        plt.axis('tight')
        plt.ylim([-5e-16,5e-16])
        #plt.xlim([(Hbeta+OIII)/2.0-0.05,(Hbeta+OIII)/2.0+0.05])
        fig = plt.gcf()
        lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2)
        fig.set_size_inches(18.5, 5.5)
        plt.setp(line1, linewidth=0.1, color='b')
        plt.setp(line2, linewidth=0.1)
        plt.setp(line3, linewidth=0.1)
        plt.setp(line5, linewidth=0.1)
        #plt.setp(line6, linewidth=0.1)
        print '../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets_specific_regions_exposures_averaged_IZ/'+hdulist_sci[i][j][-66:-57]+'_Target_4_specific_region.pdf'
        plt.savefig('../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets_specific_regions_exposures_averaged_IZ/'+hdulist_sci[i][j][-66:-57]+'_Target_4_specific_region.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')
        plt.close()
        #z = 0.869
        #if k == 20:
        #    z = 0.867
        #if k == 24:
        #    z = 0.861
        #if k == 9:
        z = 0.868
        y_range = (min(sky_subtracted_sum_data_tot2),max(sky_subtracted_sum_data_tot2))
        Hbeta = 0.4861 * ( 1 + z )
        OIII = 0.5007 * ( 1 + z )
        line1, = plt.plot(x,sky_subtracted_sum_data_tot2, label='Total Flux of Target IFU')
        line2, = plt.plot((Hbeta,Hbeta),y_range)
        line3, = plt.plot((OIII,OIII),y_range)
        line5, = plt.plot((0.78,1.08985),(0.0,0.0))
        #line6, = plt.plot(x, avg_flux_sky_list, label='Mean Flux of Scaled Sky IFU')
        plt.xlabel('Wavelength (microns)',fontsize=10)
        plt.ylabel('Flux (ergs cm$\mathregular{^{-2}}$ $\\AA\mathregular{^{-1}}$ s$\mathregular{^{-1}}$)',fontsize=10)
        plt.axis('tight')
        plt.ylim([-5e-16,5e-16])
        #plt.xlim([(Hbeta+OIII)/2.0-0.05,(Hbeta+OIII)/2.0+0.05])
        fig = plt.gcf()
        lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2)
        fig.set_size_inches(18.5, 5.5)
        plt.setp(line1, linewidth=0.1, color='b')
        plt.setp(line2, linewidth=0.1)
        plt.setp(line3, linewidth=0.1)
        plt.setp(line5, linewidth=0.1)
        #plt.setp(line6, linewidth=0.1)
        print '../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets_specific_regions_exposures_averaged_IZ/'+hdulist_sci[i][j][-66:-57]+'_Target_15_specific_region.pdf'
        plt.savefig('../figures/laptop/target_sky_fraction_of_object/sky_subtracted_spectrum_of_IFUs/specific_targets_specific_regions_exposures_averaged_IZ/'+hdulist_sci[i][j][-66:-57]+'_Target_15_specific_region.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')

        plt.close()


