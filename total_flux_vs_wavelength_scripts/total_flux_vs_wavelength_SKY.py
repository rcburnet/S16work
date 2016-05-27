#!/usr/bin/python

####
#scripts that will plot sky flux for each arm in a single sky interim fits file. Outputs are each individual arm's flux vs wavelength flux in pdf format
####

import numpy as np
import scipy as sp
import glob
import os
import fnmatch
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy.io import fits

hdulist1 = fits.open('/home/rburnet/reflex/project_data/old_reflex_files/reflex_tmp_products/kmos/kmos_sci_red_1/2016-05-16T15:04:48.737/sci_interim_sky_KMOS.2014-07-11T06:34:52.380.fits')

for i in range(len(hdulist1)):

	x = []
        total_flux = []
	
	try:
		for j in range(len(hdulist1[i].data)):
			total_flux.append(np.sum(hdulist1[i].data[j][~np.isnan(hdulist1[i].data[j])]))
			x.append(j+1)
			x[j] = 1+j*0.0001752906692721
	
		line1, = plt.plot(x,total_flux)
        #if '0034-YJ' in hdulist[j]:
                line2, = plt.plot((1.22527,1.22527),(min(total_flux),max(total_flux)))
                line3, = plt.plot((1.22927,1.22927),(min(total_flux),max(total_flux)))
                line4, = plt.plot((1.22127,1.22127),(min(total_flux),max(total_flux)))
        #else:
        #        line2, = plt.plot((1.22659,1.22659),(min(total_flux),max(total_flux)))
        #        line3, = plt.plot((1.23059,1.23059),(min(total_flux),max(total_flux)))
        #        line4, = plt.plot((1.22259,1.22259),(min(total_flux),max(total_flux)))
        	line5, = plt.plot((1.0,1.35882),(0.0,0.0))
        	line2.set_dashes([4,2,2,2])
        	line3.set_dashes([2,2])
        	line4.set_dashes([2,2])
        	plt.title('Total Flux vs Wavelength',fontsize=10)
        	plt.xlabel('Wavelength (microns)',fontsize=10)
	        plt.ylabel('Total Flux (ergs cm$\mathregular{^{-2}}$ $\\AA\mathregular{^{-1}}$ s$\mathregular{^{-1}}$)',fontsize=10)
	        plt.axis('tight')
        	fig = plt.gcf()
        	fig.set_size_inches(18.5, 5.5)
        	plt.setp(line2, linewidth=0.1, color='g')
        	plt.setp(line3, linewidth=0.1, color='r')
        	plt.setp(line4, linewidth=0.1, color='r')
        	plt.setp(line5, linewidth=0.1)
        #if not os.path.exists('./figures/'+hdulist[j][53:-48]):
        #        os.makedirs('./figures/'+hdulist[j][53:-48])
        	plt.savefig(hdulist1[i].header[63]+'.pdf')
	        plt.close()
        	hdulist[i].close()

	except:
		print "error"
