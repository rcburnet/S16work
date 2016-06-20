#!/usr/bin/python

#####
# script that will plot total flux vs wavelength while ignoring NaN of combined sky subtracted uncalibrated cubes in /home/rburnet/reflex/project_data/OBs_combined/cubes_from_raw_data/manual_sky_subtraction_cubes/ except for sky arms
#####

import numpy as np
import scipy as sp
import glob
import os
import fnmatch
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy.io import fits

hdulist=[]

for root,dirnames,filenames in os.walk('/home/rburnet/reflex/project_data/OBs_combined/cubes_from_raw_data/manual_sky_subtraction_cubes/'):
    for filename in fnmatch.filter(filenames, 'COMBINE_CUBE_OBJECT_[0-9]*'):
        hdulist.append(os.path.join(root, filename))

for j in range(len(hdulist)):

    if not os.path.isfile('./figures/manual_sky_subtracted_uncalibrated/'+hdulist[j][-46:]+'.pdf'):

        hdulist1 = fits.open(hdulist[j])
        
        if 'YJ' in hdulist[j]:
            x = np.linspace(1,1.359,2048)
        if 'IZ' in hdulist[j]:
            x = np.linspace(0.78,1.09,2048)

        total_flux = []

        for i in range(len(hdulist1[1].data)):
            total_flux.append(np.sum(hdulist1[1].data[i][~np.isnan(hdulist1[1].data[i])]))

        if not os.path.exists('./figures/manual_sky_subtracted_uncalibrated/'+hdulist[j][-46:-27]):
            os.makedirs('./figures/manual_sky_subtracted_uncalibrated/'+hdulist[j][-46:-27])
        
        line1, = plt.plot(x,total_flux)
        '''
        if '0034-YJ' in hdulist[j]:
            line2, = plt.plot((1.22527,1.22527),(min(total_flux),max(total_flux)))
            line3, = plt.plot((1.22927,1.22927),(min(total_flux),max(total_flux)))
            line4, = plt.plot((1.22127,1.22127),(min(total_flux),max(total_flux)))
        else:
            line2, = plt.plot((1.22659,1.22659),(min(total_flux),max(total_flux)))
            line3, = plt.plot((1.23059,1.23059),(min(total_flux),max(total_flux)))
            line4, = plt.plot((1.22259,1.22259),(min(total_flux),max(total_flux)))
        '''
        if 'YJ' in hdulist[j]:
            line5, = plt.plot((1.0,1.35882),(0.0,0.0))
        if 'IZ' in hdulist[j]:
            line5, = plt.plot((0.78,1.09),(0.0,0.0))
        #line2.set_dashes([4,2,2,2])
        #line3.set_dashes([2,2])
        #line4.set_dashes([2,2])
        plt.title('Total Flux vs Wavelength from \n'+hdulist[j],fontsize=10)
        plt.xlabel('Wavelength (microns)',fontsize=10)
        plt.ylabel('Total Flux (Counts)',fontsize=10)
        plt.axis('tight')
        plt.ylim([-3000,3000])
        fig = plt.gcf()
        fig.set_size_inches(18.5, 5.5)
        plt.setp(line1, linewidth=0.1, color='b')
        #plt.setp(line2, linewidth=0.1, color='g')
        #plt.setp(line3, linewidth=0.1, color='r')
        #plt.setp(line4, linewidth=0.1, color='r')
        plt.setp(line5, linewidth=0.1)
        print './figures/manual_sky_subtracted_uncalibrated/'+hdulist[j][-46:]+'.pdf'
        plt.savefig('./figures/manual_sky_subtracted_uncalibrated/'+hdulist[j][-46:]+'.pdf')
        plt.close()
        hdulist1.close()
