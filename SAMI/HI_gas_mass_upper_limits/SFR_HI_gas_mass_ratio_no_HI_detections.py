#!/usr/bin/python

# Script to calculate SFR/HI gas mass ratio for each of the 7 SAMI targets. SFR is total (integrated) SFR from the SFR datacubes created by Halpha_to_SFR.py script. HI gas mass is calculated using the first order approximation as detailed in https://iopscience.iop.org/article/10.1086/497431/pdf from the ALFALFA survey data. The calculation is not needed as the HI mass is already computed and presented in the a.70 catalog. Kept calculation in script (commented out) in case I want to look at it or need it in the future.

from astropy.io import fits
import numpy as np
import os
import fnmatch

# Import HI gas mass upper limits corresponding SAMI names of both targets without HI detections and targets with HI detections
from HI_gas_mass_upper_limits_using_50_percent_completeness_limit import HI_gas_mass_excluding_HI_detections
from HI_gas_mass_upper_limits_using_50_percent_completeness_limit import hdulist_within_ALFALFA_excluding_HI_detections
from HI_gas_mass_upper_limits_using_50_percent_completeness_limit import hdulist_within_ALFALFA_of_HI_detections

hdulist = []
filename_list = []

SFR_HI_gas_mass_ratio_list = []

for root,dirnames,filenames in os.walk('/home/rburnet/SAMI/data/SFR/with_HI_detections'):
    for filename in fnmatch.filter(filenames, '*.fits'):
        hdulist.append(os.path.join(root, filename))
        filename_list.append(filename)

for i in range(len(hdulist)):
    hdulist1 = fits.open(hdulist[i])

    data = hdulist1[0].data

    data = np.nan_to_num(data)

    HI_gas_mass = 0.0

    if '216843' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)   #Wavelength range of SAMI data
        z = 0.02382 #Redshift as detailed in SAMI-EDR catalogue
        D = 106.6   #Distance in Mpc as detailed in ALFALFA a.70 catalogue
        S = 0.99    #Peak HI flux as detailed in ALFALFA a.70 catalogue
        W = 142     #Velocity width as detailed in ALFALFA a.70 catalogue
        HI_gas_mass = 10**9.42 #9.42 from ALFALFA a.70 catalog (log(HI_gas_mass)). Units of solar mass.
    if '220371' in hdulist[i]:
        x = np.linspace(6261.12,7425.48,2048)
        z = 0.02025
        D = 90.8
        S = 1.06
        W = 187 
        HI_gas_mass = 10**9.31
    if '279917' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.01792
        D = 79.3
        S = 1.76
        W = 143
        HI_gas_mass = 10**9.42
    if '623641' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.01780
        D = 77.9
        S = 0.79
        W = 135
        HI_gas_mass = 10**9.05
    if '623726' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.01786
        D = 79.0
        S = 1.82
        W = 141
        HI_gas_mass = 10**9.43
    if '79635' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.04006
        D = 175.8
        S = 2.20
        W = 347
        HI_gas_mass = 10**10.21
    if '91924' in hdulist[i]:
        x = np.linspace(6261.12,7425.48,2048)
        z = 0.05251
        D = 228.2
        S = 1.43
        W = 63
        HI_gas_mass = 10**10.24

    theta = hdulist1[0].header['CATADEC']*np.pi/180.0   #Extract declination of source
    A = (D * 1000)**2 * (np.cos(np.pi/2.0 - theta) - np.cos(np.pi/2.0 - theta + 0.00014 * np.pi / 180.0)) * (0.00014 * np.pi / 180.0) #area in pc^2 of 1 pix. Area derived from surface element integral over 1 pix area. pi/2 - theta since declination starts from equator (pi/2), not from zenith (0.0).

    tot_SFR = np.sum(data) * A  #tot SFR is the sum of the flux elements (dSFR * A summed, or the sum of SFR * A)
    #HI_gas_mass1 = 2.356e5 * D**2.0 * S    #calculation of HI gas mass to the first order in units of solar masses.
    SFR_HI_gas_mass_ratio = tot_SFR / HI_gas_mass
    SFR_HI_gas_mass_ratio_list.append(SFR_HI_gas_mass_ratio)

SFR_HI_gas_mass_upper_limits_ratio = [] #SFR to HI gas mass ratio using upper limits of targets without HI detections

#for i in range(len(hdulist_within_ALFALFA_excluding_HI_detections)):
        

print SFR_HI_gas_mass_ratio_list
