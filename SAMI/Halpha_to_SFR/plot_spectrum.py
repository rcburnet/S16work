#Script to plot spectrum of the no HI detection 58 SAMI targets

import numpy as np
import matplotlib.pyplot as plt
from astropy.io import fits
import os
import fnmatch

from Halpha_to_SFR_and_SFR_to_HI_gas_mass_ratio_no_HI_detections import GAMA_name_list_no_HI_detections, sami_z

hdulist = []
filename_list = []

for root,dirnames,filenames in os.walk('/home/rburnet/SAMI/data/without_HI_detections/'):
    for filename in fnmatch.filter(filenames, '*.fits'):
        hdulist.append(os.path.join(root, filename))
        filename_list.append(filename)

#hdulist1 = fits.open('/home/rburnet/SAMI/data/without_HI_detections/623712_red_7_Y13SAR1_P009_09T015.fits')

FWHM_length = []    #list of lengths of FWHM_x lists for the 35 targets that have detectable Halpha lines

for i in range(len(hdulist)):
    for j in range(len(GAMA_name_list_no_HI_detections)):
        if str(GAMA_name_list_no_HI_detections[j]) in hdulist[i]:
            z = sami_z[j]
    try: 
        hdulist1 = fits.open(hdulist[i])

        centre_wavelength = hdulist1[0].header['CRVAL3']
        increment = hdulist1[0].header['CDELT3']
        x = np.linspace(centre_wavelength - increment * 1023, centre_wavelength + increment * 1024, 2048)

        halpha = 6562.8 #Angstroms

        l = halpha*(1+z) 
        NII_1 = 6548 * ( 1 + z )
        NII_2 = 6584 * ( 1 + z )
        
        total_flux = [] #Total flux of entire spectrum, to be plotted
        fwhm_flux = []  #the flux within the FWHM of the Halpha profile, a subset of total_flux

        for j in range(len(hdulist1[0].data)):
            if x[j] > NII_1 and x[j] < NII_2:
                fwhm_flux.append(np.sum(hdulist1[0].data[j][~np.isnan(hdulist1[0].data[j])]))
            total_flux.append(np.sum(hdulist1[0].data[j][~np.isnan(hdulist1[0].data[j])]))

        mean = np.mean(total_flux)

        FWHM = (max(fwhm_flux)+mean)/2.0    #The actual FWHM flux value

        FWHM_y = []
        FWHM_x = []
        for j in range(len(total_flux)):
            if total_flux[j] >= FWHM:
                if x[j] > NII_1 and x[j] < NII_2-10:    #Only place values in FWHM_y and FWHM_x between the two NII lines (-10 from NII_2 as a quick fix for 91924 source)
                    FWHM_y.append(total_flux[j])
                    FWHM_x.append(x[j])

        hrange = [] #list of True/False values. True for wavelength slices inside of FWHM range. False otherwise. To be used to choose slices to include in collapsed cube.

        #If current target is in August notes list (ie. one of the targets that has detectable Halpha lines), append the length of the FWHM_x list (the number of wavelength slices in the FWHM) to the FWHM_length list, to be used to calculated the average FWHM length at the end of the scriptto help with calculating the upper limit to the Halpha flux for the targets that don't have detectable Halpha to calculate their SFRs and thus their SFR/HI gas mass ratios for SFR_plot.py script.
        for j in ['79693', '79710', '79712', '209181', '218713', '218717', '220320', '230714', '272831', '279818', '279891', '279943', '289116', '289200', '302846', '302994', '325390', '373202', '373284', '377962', '381215', '417568', '422639', '422683', '517164', '517302', '599582', '599761', '599877', '617945', '618152', '623620', '623679', '623712', '623722']:
            if j in hdulist[i]:
                print j
                FWHM_length.append(len(FWHM_x))

        for j in range(len(x)):
            if x[j] < FWHM_x[0]:
                hrange.append(False)
            elif x[j] < FWHM_x[len(FWHM_x)-1]:
                hrange.append(True)
            else:
                hrange.append(False)

        summed_data = np.zeros([50,50])
        data = hdulist1[0].data #Doing this again to bring back the NaN values

        for j in range(0,len(data)):
            if hrange[j] == True:
                summed_data += data[j]  #summed_data is summed data within FWHM profile. ie. the collapsed Halpha image

        #Plot
        line1, = plt.plot(x,total_flux)
        line5, = plt.plot((x[0],x[2047]),(0.0,0.0))
        line6, = plt.plot((x[0],x[2047]),(mean,mean))
        line7, = plt.plot((x[0],x[2047]),(FWHM,FWHM))
        line8, = plt.plot((FWHM_x[0],FWHM_x[0]),(min(total_flux),max(total_flux)))
        line9, = plt.plot((FWHM_x[len(FWHM_x)-1],FWHM_x[len(FWHM_x)-1]),(min(total_flux),max(total_flux)))
        line2, = plt.plot((l,l),(min(total_flux),max(total_flux)))
        line3, = plt.plot((NII_1,NII_1),(min(total_flux),max(total_flux)))
        line4, = plt.plot((NII_2,NII_2),(min(total_flux),max(total_flux)))
        line2.set_dashes([4,2,2,2])
        line3.set_dashes([2,2])
        line4.set_dashes([2,2])
        plt.title('Total Flux vs Wavelength from',fontsize=10)
        plt.xlabel('Wavelength ($\\AA$)',fontsize=10)
        plt.ylabel('Total Flux (ergs cm$\mathregular{^{-2}}$ $\\AA\mathregular{^{-1}}$ s$\mathregular{^{-1}}$ pix$^{-1}$)',fontsize=10)
        plt.axis('tight')
        plt.xlim([NII_1,NII_2])
        plt.ylim([0,max(fwhm_flux)])
        fig = plt.gcf()
        fig.set_size_inches(18.5, 5.5)
        plt.setp(line1, linewidth=0.1, color='b')
        plt.setp(line2, linewidth=0.1, color='g')
        plt.setp(line3, linewidth=0.1, color='r')
        plt.setp(line4, linewidth=0.1, color='r')   
        plt.setp(line5, linewidth=0.1)
        plt.setp(line6, linewidth=0.1)
        plt.setp(line7, linewidth=0.1)
        plt.setp(line8, linewidth=0.1)
        plt.setp(line9, linewidth=0.1)

        print filename_list[i], z
        plt.savefig('figures/spectra/'+filename_list[i]+'.pdf')
        plt.close()
        hdulist1.close()
    except:
        print filename_list[i]+' didn\'t work'

#print FWHM_length, len(FWHM_length), np.mean(FWHM_length)
