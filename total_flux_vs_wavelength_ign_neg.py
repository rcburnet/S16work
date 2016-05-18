#!/usr/bin/python

#####
# script that will plot total flux vs wavelength while ignoring NaN, negative values, and outliers (> 1*10^-18), and
# plot std. dev. vs wavelength
#####

import numpy as np
import scipy as sp
import matplotlib.pyplot as plt
from astropy.io import fits

hdulist1 = fits.open('/home/rcburnet/work/project_data/reflex_end_products/2016-05-16T09:20:10/KMOS.2014-07-11T06:34:52.380_tpl/CL0034-YJ-OB-1_COMBINE_SCI_RECONSTRUCTED_11.fits')

#hdulist1 = fits.open('/home/rcburnet/work/project_data/reflex_tmp_products/kmos/kmos_sci_red_1/2016-05-16T15:04:48.737/sci_interim_sky_KMOS.2014-07-11T06:34:52.380.fits')

#hdulist1 = fits.open('/home/rcburnet/work/demo/reflex_end_products/2016-05-09T14:50:14/KMOS.2013-06-30T23:48:06.049_tpl/SCI-GUM43_COMBINE_SCI_RECONSTRUCTED_001.fits')

x = []
total_flux = []
std_flux = []

for i in range(len(hdulist1[1].data)):
	std_flux.append(np.std(hdulist1[1].data[i][~np.isnan(hdulist1[1].data[i])]))
	#### first time, ignore negatives and extreme values larger than 1e-18
	#test = hdulist1[1].data[i][hdulist1[1].data[i] >= 0]
	#test = test[test <= 1e-18]
	#### second time, ignore values that are larger than mean+#*stdev and smaller than mean-#*stddev 
	test = hdulist1[1].data[i][~np.isnan(hdulist1[1].data[i])]
	test = test[test >= np.average(test) - 2*std_flux[i]]
	test = test[test <= np.average(test) + 2*std_flux[i]]
	####
	total_flux.append(np.sum(test[~np.isnan(test)]))
	x.append(i+1)
for i in range(len(x)):
	x[i] = 1+i*0.0001752906692721

f, (ax1, ax2) = plt.subplots(1, 2, sharey=True)
ax1.plot(x,total_flux)
ax2.plot(x,std_flux)
ax1.set_title('Total Flux vs Wavelength')
ax2.set_title('Standard Dev. vs Wavelength')
ax1.set_xlabel('Wavelength (microns)')
ax2.set_xlabel('Wavelength (microns)')
ax1.set_ylabel('Total Flux for the slice')
ax2.set_ylabel('Standard Dev. for the slice')
plt.show()
