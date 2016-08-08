#!/usr/bin/python

# Script to convert Halpha to SFR of 58 SAMI sources that are within ALFALFA survey area but do not have HI detections. First, calculate mean continuum, subtract Halpha flux to continuum, integrate over Halpha profile wavelength range to get total flux for every pixel, convert that flux to luminosity using lumonisity distances (calculated using astropy.cosmology module), then use K98 Halpha luminosity to SFR relation to generate SFR map of sources.

#Script also calciulates HI gas mass upper limits of 58 SAMI targets and calculates the SFR/HI gas mass upper limits value and prints them to terminal.

# ***MUST BE RUN IN WORKING DIRECTORY WHERE YOU WANT TO PUT FITS FILES IN

from astropy.io import fits
import numpy as np
import os
import fnmatch
import matplotlib.pyplot as plt
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)   #Planck 2015 parameters


#Retrieve filenames of raw SAMI data cubes
hdulist = []
filename_list = []

for root,dirnames,filenames in os.walk('/home/rburnet/SAMI/data/without_HI_detections/'):  #change this to directory with raw fits files
    for filename in fnmatch.filter(filenames, '[0-9]*red*[0-9].fits'):
        hdulist.append(os.path.join(root, filename))
        filename_list.append(filename)


#Extract GAMA id's from filenames
GAMA_name_list = []

for i in range(len(filename_list)):
    if '_' in filename_list[i][0:6]:
        GAMA_name_list.append(filename_list[i][0:5])
    else:
        GAMA_name_list.append(filename_list[i][0:6])

#Extract ALFALFA survey a.70 velocity widths
alfalfa = open('../a70_160624.csv')
alfalfa_lines = alfalfa.readlines()

for i in range(len(alfalfa_lines)):
    alfalfa_lines[i] = alfalfa_lines[i].split(',')

alfalfa.close()

alfalfa_W = []  #ALFALFA 50% velocity widths

for i in range(len(alfalfa_lines)):
    if i != 0:
        alfalfa_W.append(float(alfalfa_lines[i][7]))

ave_W = np.mean(alfalfa_W)  #average 50% velocity width of ALFALFA survey detections

#Extract required data (coordinates and redshifts of SAMI targets)
sami = open('/home/rburnet/S16work/SAMI/SAMI_EarlyDataRelease.txt')
sami_lines = sami.readlines()

for i in range(len(sami_lines)):
    sami_lines[i] = sami_lines[i].split(' ')

sami.close()

sami_coord = []
sami_z = []

for j in range(len(GAMA_name_list)):
    for i in range(len(sami_lines)):
        if i != 0 and i != 108:
            if GAMA_name_list[j] in sami_lines[i]:
                try:
                    sami_coord.append((float(sami_lines[i][3]),float(sami_lines[i][7])))
                    try:
                        sami_z.append(float(sami_lines[i][16]))
                    except:
                        sami_z.append(float(sami_lines[i][17]))
                except:
                    sami_coord.append((float(sami_lines[i][3]),float(sami_lines[i][8])))
                    try:
                        sami_z.append(float(sami_lines[i][16]))
                    except:
                        sami_z.append(float(sami_lines[i][17]))


#Calculate luminosity distance of SAMI targets using redshifts
sami_D = []

for i in range(len(sami_z)):
    sami_D.append(cosmo.luminosity_distance(sami_z[i]).value)

x = np.linspace(6843.01442183-0.568812591597*1023,6843.01442183+0.568812591597*1024,2048)

#Create SFR images
for i in range(len(hdulist)):
    hdulist1 = fits.open(hdulist[i])
    centre_wavelength = hdulist1[0].header['CRVAL3']
    increment = hdulist1[0].header['CDELT3']
    x = np.linspace(centre_wavelength - increment * 1023, centre_wavelength + increment * 1024, 2048)
    D = sami_D[i] * 3.0856776e+24  #Mpc to cm
    z = sami_z[i]

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

    #Now calculate SFR to HI gas mass ratios using upper limits to HI gas mass

    #Calculate HI gas mass upper limits of SAMI targets
    log_S_21_90 = 0.5 * np.log10(ave_W) - 1.14   #log of 90% completeness limit integrated flux density from Haynes et al 2011 (https://arxiv.org/pdf/1109.0027v1.pdf)

    log_S_21_50 = log_S_21_90 - 0.067   #log of 50% completeness limit integrated flux density

    S_21_50 = 10**(log_S_21_50) #50% completeness limit integrated flux density, units are in Jy km/s

    HI_gas_mass = 2.356e5 * (D / 3.0856776e+24)**2.0 * S_21_50  #upper limit to HI gas mass of SAMI targets using Giovanelli et al, 2005 relations. Units are in solar masses
    data = SFR

    data = np.nan_to_num(data)

    tot_SFR = np.sum(data) * A  #tot SFR is the sum of the flux elements (dSFR * A summed, or the sum of SFR * A)
    SFR_HI_gas_mass_ratio = tot_SFR / HI_gas_mass
    print hdulist[i], SFR_HI_gas_mass_ratio #prints name of fits file processed and the SFR/HI gas mass upper limits ratio (corresponds to lower limit on ratio)

    hdulist1.close()

