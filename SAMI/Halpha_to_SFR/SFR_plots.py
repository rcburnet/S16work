#Script that creates two plots - SFR/HI gas mass vs Mstar plot and SFR vs HI gas mass plot - of the 65 SAMI-EDR targets in the ALFALFA survey area (7 with HI detections, and thus calculable HI gas mass and SFR/HI gas mass ratios, and 58 with no HI detections, and thus no calculable HI gas mass except for HI gas mass upper limits using 50% survey completeness limit which is used for the SFR/HI gas mass ratios).

import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams['xtick.labelsize'] = 10    #xtick label size
mpl.rcParams['ytick.labelsize'] = 10    #ytick label size

#import GAMA name lists, SFR/HI gas mass ratios lists, HI gas mass lists, SFR lists, and Mstar lists of 65 targets. NOTE: Mstar is stellar mass within SAMI aperture.
print '7 targets with HI detection \n'
from SFR_HI_gas_mass_ratio import GAMA_name_list_with_HI_detections, SFR_HI_gas_mass_ratio_list_with_HI_detections, HI_gas_mass_list_with_HI_detections, tot_SFR_list_with_HI_detections, sami_Mstar_list_with_HI_detections

print '\n'
print '58 targets with no HI detection \n'
from Halpha_to_SFR_and_SFR_to_HI_gas_mass_ratio_no_HI_detections import GAMA_name_list_no_HI_detections, SFR_HI_gas_mass_ratio_list_no_HI_detections_detectable_Halpha, SFR_HI_gas_mass_ratio_list_no_HI_detections_no_detectable_Halpha, HI_gas_mass_list_no_HI_detections_detectable_Halpha, tot_SFR_list_no_HI_detections_detectable_Halpha, HI_gas_mass_list_no_HI_detections_no_detectable_Halpha, tot_SFR_list_no_HI_detections_no_detectable_Halpha, Mstar_list_no_HI_detections_detectable_Halpha, Mstar_list_no_HI_detections_no_detectable_Halpha

slope, intercept, r_value, p_value, std_err = stats.linregress(np.log10(HI_gas_mass_list_no_HI_detections_detectable_Halpha), np.log10(tot_SFR_list_no_HI_detections_detectable_Halpha))
r_squared = r_value**2.0
best_fit = slope * np.log10(HI_gas_mass_list_no_HI_detections_detectable_Halpha) + intercept

print slope, intercept

#Plot plots

#plot SFR vs HI gas mass
fig = plt.figure()
plt.scatter(np.log10(HI_gas_mass_list_with_HI_detections), np.log10(tot_SFR_list_with_HI_detections), label = 'Targets with HI detections', color = 'red')
plt.scatter(np.log10(HI_gas_mass_list_no_HI_detections_detectable_Halpha), np.log10(tot_SFR_list_no_HI_detections_detectable_Halpha), label = 'Targets with no HI detections, \ndetectable Halpha', color = 'blue')
plt.scatter(np.log10(HI_gas_mass_list_no_HI_detections_no_detectable_Halpha), np.log10(tot_SFR_list_no_HI_detections_no_detectable_Halpha), label = 'Targets with no HI detections, \nno detectable Halpha', color = 'green')
#plt.plot(np.log10(HI_gas_mass_list_no_HI_detections_detectable_Halpha), best_fit, color = 'black')
plt.title('SFR vs HI gas mass of 65 SAMI targets in ALFALFA survey area', fontsize=10)
plt.ylabel('Log SFR (M$_\odot$/yr)', fontsize=10)
plt.xlabel('Log M$_{gas}$ (M$_\odot$)', fontsize=10)
plt.axis('tight')
lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2, fontsize=10)
plt.savefig('SFR_vs_HI_gas_mass.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')
plt.close()

#plot SFR/HI gas mass ratio vs Mstar
fig = plt.figure()
plt.scatter(np.log10(sami_Mstar_list_with_HI_detections), np.log10(SFR_HI_gas_mass_ratio_list_with_HI_detections), label = 'Targets with HI detections', color = 'red')
plt.scatter(np.log10(Mstar_list_no_HI_detections_detectable_Halpha), np.log10(SFR_HI_gas_mass_ratio_list_no_HI_detections_detectable_Halpha), label = 'Targets with no HI detections, \ndetectable Halpha', color = 'blue')
plt.scatter(np.log10(Mstar_list_no_HI_detections_no_detectable_Halpha), np.log10(SFR_HI_gas_mass_ratio_list_no_HI_detections_no_detectable_Halpha), label = 'Targets with no HI detections, \nno detectable Halpha', color = 'green')
plt.title('SFR/HI gas mass ratio vs stellar mass of 65 SAMI targets in ALFALFA survey area', fontsize=10)
plt.ylabel('Log SFR/M$_{gas}$ (yr$^{-1}$)', fontsize=10)
plt.xlabel('Log M$_{*}$ (M$_\odot$)', fontsize=10)
plt.axis('tight')
lgd = plt.legend(bbox_to_anchor=(1.02, 1.02), loc=2, fontsize=10)
plt.savefig('SFR_HI_gas_mass_ratio_vs_Mstar.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')
plt.close()

