#!/usr/bin/python

# Script to convert Halpha to SFR of 7 SAMI sources. First, calculate mean continuum, subtract Halpha flux to continuum, integrate over Halpha profile wavelength range to get total flux for every pixel, convert that flux to luminosity with known distances from ALFALFA csv files, then use K98 Halpha luminosity to SFR relation to generate SFR map of sources.

# ***MUST BE RUN IN WORKING DIRECTORY WHERE YOU WANT TO PUT FITS FILES IN

from astropy.io import fits
import numpy as np
import os
import fnmatch
import matplotlib.pyplot as plt

hdulist = []
filename_list = []

for root,dirnames,filenames in os.walk('/home/rburnet/SAMI/data/with_HI_detections/'):  #change this to directory with raw fits files
    for filename in fnmatch.filter(filenames, '[0-9]*red*[0-9].fits'):
        hdulist.append(os.path.join(root, filename))
        filename_list.append(filename)

for i in range(len(hdulist)):
    hdulist1 = fits.open(hdulist[i])

    if '216843' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.02382
        D = 106.6 * 3.0856776e+24  #Mpc to cm
    if '220371' in hdulist[i]:
        x = np.linspace(6261.12,7425.48,2048)
        z = 0.02025
        D = 90.8 * 3.0856776e+24
    if '279917' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.01792
        D = 79.3 * 3.0856776e+24
    if '623641' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.01780
        D = 77.9 * 3.0856776e+24
    if '623726' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.01786
        D = 79.0 * 3.0856776e+24
    if '79635' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.04006
        D = 175.8 * 3.0856776e+24
    if '91924' in hdulist[i]:
        x = np.linspace(6261.12,7425.48,2048)
        z = 0.05251
        D = 228.2 * 3.0856776e+24

    halpha = 6562.8 #Angstroms

    l = halpha*(1+z)
    NII_1 = 6548 * ( 1 + z )
    NII_2 = 6584 * ( 1 + z )

    Halpha_flux = []    #the flux within the Halpha profile
    continuum_flux = [] #the continuum flux (between Halpha+7 and NII_2-7)
    
    # retrieve continuum flux
    data = hdulist1[0].data * 1e-16
    for j in range(len(data)):
        if x[j] > l+7 and x[j] < NII_2-7:    #Only place values between Halpha+7 and NII_2-7 (to get continuum)
            continuum_flux.append(data[j])

    #continuum is mean over continuum wavelength range for each pixel, 50x50 array where each pixel is the continuum flux mean for that specific pixel.
    continuum_flux = np.mean(continuum_flux, axis=0)
 
    # subtract Halpha profile from continuum, then integrate over Halpha profile to get Halpha flux
    for j in range(len(data)):
        if x[j] > l-7 and x[j] < l+7:   #look over Halpha profile and a little in the continuum to subtract continuum from flux and then integrate
            Halpha_flux.append(data[j] - continuum_flux)

    Halpha_flux = np.array(Halpha_flux)

    #loop into Halpha_flux array to integrate to get true Halpha flux for every pixel. Use trapezoid method to integrate.
    dx = x[1]-x[0]
    integral = np.zeros([50,50])
    for j in range(len(Halpha_flux)):
        if j != 0 and j != len(Halpha_flux)-1:
            Halpha_flux[j] = Halpha_flux[j]*2.0 #All f(x) between first and last point must be doubled, see wikipedia article on trapezoid rule on why.
        integral += Halpha_flux[j]

    integral * dx / 2.0  #final integrated array

    #convert Halpha_flux to luminosity
    L = 4 * np.pi * D**2.0 * integral

    #Convert luminosity to SFR
    theta = hdulist1[0].header['CATADEC']*np.pi/180.0   #Extract declination of source
    A = (D * 1000 / 3.0856776e+24)**2 * (np.cos(np.pi/2.0 - theta) - np.cos(np.pi/2.0 - theta + 0.00014 * np.pi / 180.0)) * (0.00014 * np.pi / 180.0) #area in pc^2 of 1 pix. Area derived from surface element integral over 1 pix area. pi/2 - theta since declination starts from equator (pi/2), not from zenith (0.0).
    SFR = 7.9e-42 * L  / A #Calculate SFR from luminosity, convert from SFR/pix to SFR/pc^2.

    hdulist1[0].data = SFR
    hdulist1[0].header['BUNIT'] = 'M_sun /yr /pc^2'
    hdulist1.writeto(filename_list[i]+'SFR.fits')
    hdulist1.close()

