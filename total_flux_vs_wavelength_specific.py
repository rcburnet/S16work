#!/usr/bin/python

#####
# script that will plot total flux vs wavelength while ignoring NaN for particularly bright targets (targets we expect to see high H-alpha emission)
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

for root,dirnames,filenames in os.walk('/home/rburnet/reflex/project_data/reflex_end_products/'):
	for filename in fnmatch.filter(filenames, 'CL0034-YJ-OB-[1|2]_COMBINE_SCI_RECONSTRUCTED_5.fits'):
		hdulist.append(os.path.join(root, filename))
        for filename in fnmatch.filter(filenames, 'CL0034-YJ-OB-[1|2]_COMBINE_SCI_RECONSTRUCTED_11.fits'):
                hdulist.append(os.path.join(root, filename))
        for filename in fnmatch.filter(filenames, 'CL0036-YJ-OB-[1|2]_COMBINE_SCI_RECONSTRUCTED_4.fits'):
                hdulist.append(os.path.join(root, filename))
        for filename in fnmatch.filter(filenames, 'CL0036-YJ-OB-[1|2]_COMBINE_SCI_RECONSTRUCTED_15.fits'):
                hdulist.append(os.path.join(root, filename))

for j in range(len(hdulist)):

	if not os.path.isfile('./figures/laptop/specific_targets/'+hdulist[j][53:]+'.pdf'):

		hdulist1 = fits.open(hdulist[j])

		x = []
		total_flux = []
		std_flux = []

		for i in range(len(hdulist1[1].data)):
			total_flux.append(np.sum(hdulist1[1].data[i][~np.isnan(hdulist1[1].data[i])]))
			std_flux.append(np.std(hdulist1[1].data[i][~np.isnan(hdulist1[1].data[i])]))
			x.append(i+1)
        		x[i] = 1+i*0.0001752906692721

		if not os.path.exists('./figures/laptop/specific_targets/'+hdulist[j][53:-48]):
        	        os.makedirs('./figures/laptop/specific_targets/'+hdulist[j][53:-48])
		
		line1, = plt.plot(x,total_flux)
		if 'RECONSTRUCTED_5' in hdulist[j]:
                        z = 0.869
                if 'RECONSTRUCTED_11' in hdulist[j]:
                        z = 0.867
                if 'RECONSTRUCTED_4' in hdulist[j]:
                        z = 0.861
                if 'RECONSTRUCTED_15' in hdulist[j]:
                        z = 0.868
                l = 0.6563 * ( 1 + z )
                NII_1 = 0.6548 * ( 1 + z )
                NII_2 = 0.6584 * ( 1 + z )
                line2, = plt.plot((l,l),(min(total_flux),max(total_flux)))
                line3, = plt.plot((NII_1,NII_1),(min(total_flux),max(total_flux)))
                line4, = plt.plot((NII_2,NII_2),(min(total_flux),max(total_flux)))
                line5, = plt.plot((1.0,1.35882),(0.0,0.0))
                line2.set_dashes([4,2,2,2])
                line3.set_dashes([2,2])
                line4.set_dashes([2,2])
                plt.title('Total Flux vs Wavelength from \n'+hdulist[j],fontsize=10)
		plt.xlabel('Wavelength (microns)',fontsize=10)
		plt.ylabel('Total Flux (ergs cm$\mathregular{^{-2}}$ $\\AA\mathregular{^{-1}}$ s$\mathregular{^{-1}}$)',fontsize=10)
		plt.axis('tight')
                plt.xlim([l-0.05,l+0.05])
		fig = plt.gcf()
		fig.set_size_inches(18.5, 5.5)
		plt.setp(line1, linewidth=0.1, color='b')
		plt.setp(line2, linewidth=0.1, color='g')
		plt.setp(line5, linewidth=0.1)
		print './figures/laptop/specific_targets/'+hdulist[j][53:]+'.pdf'
		plt.savefig('./figures/laptop/specific_targets/'+hdulist[j][53:]+'.pdf')
		plt.close()
		hdulist1.close()
