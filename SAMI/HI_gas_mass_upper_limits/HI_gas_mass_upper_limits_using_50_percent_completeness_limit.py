# Script to calculate upper limit to HI gas mass of the 58 SAMI targets in the ALFALFA survey area that don't have (yet or ever from ALFALFA) HI detections associated with them in the survey data. Calculates the 50% completeness limit from code 1 sources with equations 4 and 5 found in Haynes et al 2011 (https://arxiv.org/pdf/1109.0027v1.pdf) and then uses that flux density to calculate an upper limit to each of the target's HI gas mass using the HI gas mass formula found in Giovanelli et al, 2005 (https://arxiv.org/pdf/astro-ph/0508301v1.pdf). Prints results to terminal and outputs to text file "HI_gas_mass_upper_limits.txt".

#Note: This prints and outputs expected upper limit of HI mass for ALL SAMI targets in ALFALFA survey area (including those that ALFALFA has detected HI for), therefore must ignore those targets which have detected HI in the ALFALFA survey as we don't need to calculate the upper limit to the HI gas mass for those targets (as we already have their actual HI gas mass). I've kept them in for diagnostic purposes, to see if the actual HI gas masses are greater than the upper limits calculated here (expected). Of course, this is assuming they all have the average a.70 survey 50% velocity width and rms noise values; some targets will have less than the average values, resulting in perhaps less actual HI gas mass than calculated upper limit which assumes average a.70 values.
 

import numpy as np
from astropy.cosmology import FlatLambdaCDM
cosmo = FlatLambdaCDM(H0=67.8, Om0=0.308)   #Planck 2015 parameters


#Extract required data (coordinates and redshifts of SAMI targets and average 50% velocity width of ALFALFA detections)
sami = open('../SAMI_EarlyDataRelease_modified.txt')
sami_lines = sami.readlines()

for i in range(len(sami_lines)):
    sami_lines[i] = sami_lines[i].split(' ')

sami.close()

alfalfa = open('../a70_160624.csv')
alfalfa_lines = alfalfa.readlines()

for i in range(len(alfalfa_lines)):
    alfalfa_lines[i] = alfalfa_lines[i].split(',')

alfalfa.close()

sami_coord = [] #SAMI coordinates

hdulist = []    #SAMI names

for i in range(len(sami_lines)):
    if i != 0 and i != 108:
        sami_coord.append((float(sami_lines[i][1]),float(sami_lines[i][2])))
        hdulist.append(sami_lines[i][0])

alfalfa_W = []  #ALFALFA 50% velocity widths
alfalfa_coord = []

for i in range(len(alfalfa_lines)):
    if i != 0:
        alfalfa_coord.append((float(alfalfa_lines[i][4]),float(alfalfa_lines[i][5])))
        alfalfa_W.append(float(alfalfa_lines[i][7]))

ave_W = np.mean(alfalfa_W)  #average 50% velocity width of ALFALFA survey detections

#Find SAMI targets that are within ALFALFA survey area
new_closest = []    #Array to hold closest pairs between SAMI targets and ALFALFA survey data detections
within_ALFALFA = [] #Array to hold coordinates of SAMI targets that are within ALFALFA survey area, regardless of whether or not they have HI detections in the ALFALFA data
sami_z = [] #SAMI redshifts
hdulist_within_ALFALFA = [] #SAMI names of SAMI targets withing ALFALFA survey area

gama_name_list = []

for i in range(len(sami_coord)):
    new_closest.append([])
    # Find SAMI targets that are within ALFALFA range
    if (sami_coord[i][0] < 247.5 and sami_coord[i][0] > 112.5) or sami_coord[i][0] < 45.0  or sami_coord[i][0] > 330.0:
        if sami_coord[i][1] < 36.0 and sami_coord[i][1] > 0.0:
            within_ALFALFA.append(sami_coord[i])
            hdulist_within_ALFALFA.append(hdulist[i])
            sami_z.append(float(sami_lines[i+1][6]))
            gama_name_list.append(sami_lines[i+1][17])

#Calculate luminosity distance of SAMI targets using redshifts

sami_D = []

for i in range(len(sami_z)):
    sami_D.append(cosmo.luminosity_distance(sami_z[i]).value)

#Calculate HI gas mass upper limits of SAMI targets

log_S_21_90 = 0.5 * np.log10(ave_W) - 1.14   #log of 90% completeness limit integrated flux density from Haynes et al 2011 (https://arxiv.org/pdf/1109.0027v1.pdf)

log_S_21_50 = log_S_21_90 - 0.067   #log of 50% completeness limit integrated flux density

S_21_50 = 10**(log_S_21_50) #50% completeness limit integrated flux density, units are in Jy km/s

print ave_W, S_21_50

HI_gas_mass = []

#Write to text file

text_file = open('HI_gas_mass_upper_limits_using_50_percent_completeness_limit.txt', 'w')

text_file.write('GAMA name, SAMI coord, calculated HI gas mass upper limits (1e6 M_sun), actual HI gas mass (1e6 M_sun) \n')

hdulist_within_ALFALFA_of_HI_detections = []    #list of SAMI names of targets with known HI detections

for i in range(len(sami_D)):
    HI_gas_mass.append(2.356e5 * sami_D[i]**2.0 * S_21_50/1e6)  #upper limit to HI gas mass of SAMI targets using Giovanelli et al, 2005 relations, units are 1e6 M_sun
    print within_ALFALFA[i], '%i' % HI_gas_mass[i] #print coordinates of SAMI targets in ALFALFA survey area and associated upper limit HI gas mass
    if 139.77630112 == within_ALFALFA[i][0]:
        text_file.write(str(gama_name_list[i])+', '+str(within_ALFALFA[i]) + ', ' + str('%i' % HI_gas_mass[i]) + ', ' + '1120' + '\n')
        hdulist_within_ALFALFA_of_HI_detections.append(hdulist_within_ALFALFA[i])
    elif 139.99533336 == within_ALFALFA[i][0]:
        text_file.write(str(gama_name_list[i])+', '+str(within_ALFALFA[i]) + ', ' + str('%i' % HI_gas_mass[i]) + ', ' + '2630' + '\n')
        hdulist_within_ALFALFA_of_HI_detections.append(hdulist_within_ALFALFA[i])
    elif 140.09157464 == within_ALFALFA[i][0]:
        text_file.write(str(gama_name_list[i])+', '+str(within_ALFALFA[i]) + ', ' + str('%i' % HI_gas_mass[i]) + ', ' + '2690' + '\n')
        hdulist_within_ALFALFA_of_HI_detections.append(hdulist_within_ALFALFA[i])
    elif 140.19242092 == within_ALFALFA[i][0]:
        text_file.write(str(gama_name_list[i])+', '+str(within_ALFALFA[i]) + ', ' + str('%i' % HI_gas_mass[i]) + ', ' + '2630' + '\n')
        hdulist_within_ALFALFA_of_HI_detections.append(hdulist_within_ALFALFA[i])
    elif 181.237152 == within_ALFALFA[i][0]:
        text_file.write(str(gama_name_list[i])+', '+str(within_ALFALFA[i]) + ', ' + str('%i' % HI_gas_mass[i]) + ', ' + '2040' + '\n')
        hdulist_within_ALFALFA_of_HI_detections.append(hdulist_within_ALFALFA[i])
    elif 214.25699075 == within_ALFALFA[i][0]:
        text_file.write(str(gama_name_list[i])+', '+str(within_ALFALFA[i]) + ', ' + str('%i' % HI_gas_mass[i]) + ', ' + '17380' + '\n')
        hdulist_within_ALFALFA_of_HI_detections.append(hdulist_within_ALFALFA[i])
    elif 222.51551376 == within_ALFALFA[i][0]:
        text_file.write(str(gama_name_list[i])+', '+str(within_ALFALFA[i]) + ', ' + str('%i' % HI_gas_mass[i]) + ', ' + '16220' + '\n')
        hdulist_within_ALFALFA_of_HI_detections.append(hdulist_within_ALFALFA[i])
    else:
        text_file.write(str(gama_name_list[i])+', '+str(within_ALFALFA[i]) + ', ' + str('%i' % HI_gas_mass[i]) + '\n')

hdulist_within_ALFALFA_excluding_HI_detections = [] #list of SAMI names excluding the targets with known HI detections
HI_gas_mass_excluding_HI_detections = []    #HI gas mass upper limits excluding the targets with known HI detections

for i in range(len(sami_D)):
    if 139.77630112 != within_ALFALFA[i][0] and 139.99533336 != within_ALFALFA[i][0] and 140.09157464 != within_ALFALFA[i][0] and 140.19242092 != within_ALFALFA[i][0] and 181.237152 != within_ALFALFA[i][0] and 214.25699075 != within_ALFALFA[i][0] and 222.51551376 != within_ALFALFA[i][0]:
        HI_gas_mass_excluding_HI_detections.append(2.356e5 * sami_D[i]**2.0 * S_21_50/1e6)
        hdulist_within_ALFALFA_excluding_HI_detections.append(hdulist_within_ALFALFA[i])

text_file.close()

#Note: This prints and outputs expected upper limit of HI flux for ALL SAMI targets in ALFALFA survey area (including those that ALFALFA has detected HI for), therefore must ignore those targets which have detected HI in the ALFALFA survey as we don't need to calculate the upper limit to the HI gas mass for those targets (as we already have their actual HI gas mass). I've kept them in for diagnostic purposes, to see if the actual HI gas masses are greater than the upper limits calculated here (expected). Of course, this is assuming they all have the average 50% velocity width and rms noise values; some targets have less, resulting in perhaps less actual Hi gas mass than calculated upper limit.
