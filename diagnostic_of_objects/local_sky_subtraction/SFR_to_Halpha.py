#!/usr/bin/python

# Script that will take SFR(OII) values from Seans catalogue (txt tables) and convert it to Halpha flux using K98 relation (SFR = 7.9*10e-42*L(Halpha), L(Halpha) = 4*pi*D**2*Flux(Halpha)). Does it for both corrected SFR(OII) values and uncorrected values.

import numpy as np
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)   #Planck 2015 parameters


### Read txt table files for target/arm info
file_CL0034 = open('/home/rburnet/S16work/txt_tables/KMOS_GCLASS_CLUS0034.txt','r')
file_CL0036 = open('/home/rburnet/S16work/txt_tables/KMOS_GCLASS_CLUS0036_YJ.txt','r')

CL0034_info = file_CL0034.readlines()
CL0036_info = file_CL0036.readlines()

file_CL0034.close()
file_CL0036.close()

for i in range(len(CL0034_info)):
    CL0034_info[i] = CL0034_info[i].split(' ')
    for j in range(len(CL0034_info[i])):
        if j < 3:
            CL0034_info[i][j] = int(CL0034_info[i][j])
    else:
            CL0034_info[i][j] = float(CL0034_info[i][j])

for i in range(len(CL0036_info)):
    CL0036_info[i] = CL0036_info[i].split(' ')
    for j in range(len(CL0036_info[i])):
        if j < 3:
            CL0036_info[i][j] = int(CL0036_info[i][j])
        else:
            CL0036_info[i][j] = float(CL0036_info[i][j])


#Define lists for specific values needed for conversions
CL0034_redshift_values = []
CL0034_SFR_values = []
CL0034_SFR_uncor_values = []
CL0034_D = []
CL0036_redshift_values = []
CL0036_SFR_values = []
CL0036_SFR_uncor_values = []
CL0036_D = []


#Extracts needed values from txt files and determine luminosity distances
for i in range(len(CL0034_info)):
    CL0034_redshift_values.append(float(CL0034_info[i][5]))                     #Redshift values, column 6 in txt files
    CL0034_SFR_values.append(float(CL0034_info[i][17]))                         #SFR values, column 18 in txt files
    CL0034_SFR_uncor_values.append(float(CL0034_info[i][18]))                   #SFR_uncorrected values, column 19 in txt files
    CL0034_D.append(cosmo.luminosity_distance(CL0034_redshift_values[i]).value) #Luminosity distances calculated using astropy.cosmology

for i in range(len(CL0036_info)):
    CL0036_redshift_values.append(float(CL0036_info[i][5]))
    CL0036_SFR_values.append(float(CL0036_info[i][17]))
    CL0036_SFR_uncor_values.append(float(CL0036_info[i][18]))
    CL0036_D.append(cosmo.luminosity_distance(CL0036_redshift_values[i]).value)


#Turn lists into numpy arrays for mathematical manipulations without the use of for loops
CL0034_redshift_values = np.array(CL0034_redshift_values)
CL0034_SFR_values = np.array(CL0034_SFR_values)
CL0034_SFR_uncor_values = np.array(CL0034_SFR_uncor_values)
CL0034_D = np.array(CL0034_D)
CL0036_redshift_values = np.array(CL0036_redshift_values)
CL0036_SFR_values = np.array(CL0036_SFR_values)
CL0036_SFR_uncor_values = np.array(CL0036_SFR_uncor_values)
CL0036_D = np.array(CL0036_D)


#Do conversion (SFR -> Luminosity -> total Halpha flux)
CL0034_L = CL0034_SFR_values / 7.9e-42
CL0034_Halpha = CL0034_L / (4 * np.pi * (CL0034_D * 3.0856776e+24)**2.0)        #3.0856776e+24 is conversion factor from Mpc to cm
CL0034_L_uncor = CL0034_SFR_uncor_values / 7.9e-42
CL0034_Halpha_uncor = CL0034_L_uncor / (4 * np.pi * (CL0034_D * 3.0856776e+24)**2.0)

CL0036_L = CL0036_SFR_values / 7.9e-42
CL0036_Halpha = CL0036_L / (4 * np.pi * (CL0036_D * 3.0856776e+24)**2.0)
CL0036_L_uncor = CL0036_SFR_uncor_values / 7.9e-42
CL0036_Halpha_uncor = CL0036_L_uncor / (4 * np.pi * (CL0036_D * 3.0856776e+24)**2.0)


#Print results to terminal
for i in range(len(CL0034_Halpha)):
    print 'CL0034_'+str(CL0034_info[i][0]), CL0034_Halpha[i], CL0034_Halpha_uncor[i]

for i in range(len(CL0036_Halpha)):
    print 'CL0036_'+str(CL0036_info[i][0]), CL0036_Halpha[i], CL0036_Halpha_uncor[i]
