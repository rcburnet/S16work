#!/usr/bin/python

####
#script that will plot a sky flux and object flux relative residuals on the same plot to compare the two curves to see if they properly overlay and show good results. Output is sky_flux.pdf. Inputs are interim sky frames and reduced data cube of object without sky subtraction. Very messy, needs cleaning.
####

import numpy as np
import scipy as sp
import glob
import os
import fnmatch
import matplotlib.pyplot as plt
from matplotlib import rc
from astropy.io import fits

hdulist=[]	#list to hold all interim sky fits files

hdufits_of_obj_4 = []

for root,dirnames,filenames in os.walk('/home/rcburnet/work/project_data/old_reflex_files/reflex_tmp_products/kmos/kmos_sci_red_1/2016-05-16T15:04:48.737/'):
        for filename in fnmatch.filter(filenames, 'sci_interim_sky_*'):
                hdulist.append(os.path.join(root, filename))
	for filename in fnmatch.filter(filenames, 'sci_interim_object_*'):
		hdufits_of_obj_4.append(os.path.join(root, filename))

#hdufits_of_obj_4 = fits.open('/home/rcburnet/work/project_data/reflex_tmp_products/kmos/kmos_sci_red_1/2016-05-16T15:04:48.737/sci_interim_object_KMOS.2014-07-11T06:34:52.380.fits') 		##use if you want to plot idividual IFU flux

#hdufits_of_obj_4 = fits.open('/home/rcburnet/work/project_data/reflex_end_products/sky_tweak_FALSE_no_subtract_TRUE/2016-05-13T15:03:41/KMOS.2014-07-11T06:34:52.380_tpl/CL0034-YJ-OB-1_COMBINE_SCI_RECONSTRUCTED_4.fits')	##use if you want to plot specific target after data reduction without sky subtraction

total_flux_tot=[]	#array of total flux values

for j in range(len(hdulist)):

        hdulist1 = fits.open(hdulist[j])	#open each fits file in hdulist list

	x = []					#wavelength array
        total_flux = []				#temporary total flux array, used for each fits file

        for i in range(len(hdulist1[20].data)):	#look at IFU1 data and loop through data array
		total_flux.append(np.sum(hdulist1[20].data[i][~np.isnan(hdulist1[20].data[i])]))	#sum all flux and append to temporary total flux array
		x.append(i+1)			#add 1 for every loop iteration to x
	
#		if i == 947:	#location of largest peak
#			print total_flux[i]
	
	for i in range(len(x)):			#loop to correct units of x
                x[i] = 1+i*0.0001752906692721
	total_flux_tot.append(total_flux)	#append temporary total flux arrays to total flux array

total_flux = []					#reset temporary total flux array
y=[]
for j in range(len(hdufits_of_obj_4)):

	hdufits_obj = fits.open(hdufits_of_obj_4[j])
	
	total_flux = []
	
	for i in range(len(hdufits_obj[20].data)):	#loop into object fits file data
		#total_flux.append(np.sum(hdufits_obj[1].data[i][~np.isnan(hdufits_obj[1].data[i])][abs(hdufits_obj[1].data[i][~np.isnan(hdufits_obj[1].data[i])])<1000]))	#sum all flux and append to temporary total flux array, only accepts values between specified threshold
		total_flux.append(np.sum(hdufits_obj[20].data[i][~np.isnan(hdufits_obj[20].data[i])]))      #sum all flux and append to temporary total flux array
		if i == 947:		#location of largest peak
			y.append(hdufits_obj[20].data[i])	
	total_flux_tot.append(total_flux)	#append temproary total flux array to total flux array

mean = []	#array to hold mean sky flux for each wavelength

for i in range(len(total_flux_tot)):	#loop into total flux array
	for j in range(len(total_flux_tot[i])):	#loop into each array in total flux array
		if total_flux_tot[i][j]!=total_flux_tot[i][j]:	#if value of array element is NaN, set it to 0
			total_flux_tot[i][j]=0

for i in range(len(total_flux_tot[0])):	#loop into each element of sky flux arrays
	mean.append(np.mean(np.array([total_flux_tot[0][i],total_flux_tot[1][i],total_flux_tot[2][i]]))) #append mean of each sky flux to mean array

for i in range(len(total_flux_tot)):	#loop into flux array
        for j in range(len(total_flux_tot[i])):	#loop into total flux array to calculate relative residuals
		if mean[j] != 0:
	                total_flux_tot[i][j] = (total_flux_tot[i][j] - mean[j])/mean[j]

#Plot
#line6, = plt.plot(x,total_flux_tot[1],label='Sky 2014-07-11T06:50:52.037')
line8, = plt.plot(x,total_flux_tot[3],label='Obj 2014-07-11T06:34:52.380')
line1, = plt.plot(x,total_flux_tot[0], label='Sky 2014-07-11T06:42:39.007')
line9 = plt.plot(x,total_flux_tot[4],label='Obj 2014-07-11T06:50:52.037')
line10, = plt.plot(x,total_flux_tot[5],label='Obj 2014-07-11T06:59:01.745')
line7, = plt.plot(x,total_flux_tot[2],label='Sky 2014-07-11T07:06:49.051')
plt.title('Total Flux vs Wavelength of Skies and Object of Target 11 (IFU 20)',fontsize=10)
plt.xlabel('Wavelength (microns)',fontsize=10)
plt.ylabel('Total Flux',fontsize=10)
plt.axis('tight')
fig = plt.gcf()
fig.set_size_inches(18.5, 5.5)
#plt.ylim(-0.5,0.5)
line1.set_dashes([1])
#line6.set_dashes([4,2,2,2])
line7.set_dashes([1])
plt.setp(line1, linewidth=0.1, color='orange')
#plt.setp(line6, linewidth=0.1, color='r')
plt.setp(line7, linewidth=0.1, color='b')
plt.setp(line8, linewidth=0.1, color='r')
plt.setp(line9, linewidth=0.1, color='y')
plt.setp(line10, linewidth=0.1, color='g')
lgd = plt.legend(bbox_to_anchor=(1.02,1.02),loc=2,numpoints=1)
plt.savefig('sky_flux.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')
hdulist1.close()
