#Script to create collapsed cube between FWHM of Halpha profile (ie. create Halpha image) and plot spectrum

from astropy.io import fits
import numpy as np
import os
import fnmatch
import matplotlib.pyplot as plt


hdulist = []
filename_list = []

for root,dirnames,filenames in os.walk('/home/rburnet/SAMI/data/'):
    for filename in fnmatch.filter(filenames, '[0-9]*red*[0-9].fits'):
        hdulist.append(os.path.join(root, filename))
        filename_list.append(filename)

for i in range(len(hdulist)):
    hdulist1 = fits.open(hdulist[i])

    data = hdulist1[0].data

    data = np.nan_to_num(data)

    if '216843' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.02382
    if '220371' in hdulist[i]:
        x = np.linspace(6261.12,7425.48,2048)
        z = 0.02025
    if '279917' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.01792
    if '623641' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.01780
    if '623726' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.01786
    if '79635' in hdulist[i]:
        x = np.linspace(6210.35,7375.84,2048)
        z = 0.04006
    if '91924' in hdulist[i]:
        x = np.linspace(6261.12,7425.48,2048)
        z = 0.05251

    halpha = 6562.8 #Angstroms

    l = halpha*(1+z)
    NII_1 = 6548 * ( 1 + z )
    NII_2 = 6584 * ( 1 + z )

    total_flux = [] #Total flux of entire spectrum, to be plotted
    fwhm_flux = []  #the flux within the FWHM of the Halpha profile, a subset of total_flux

    for j in range(len(hdulist1[0].data)):
        if x[j] > NII_1 and x[j] < NII_2:
            fwhm_flux.append(np.sum(hdulist1[0].data[j][~np.isnan(hdulist1[0].data[j])]))    
        total_flux.append(np.sum(hdulist1[0].data[j][~np.isnan(hdulist1[0].data[j])]))

    mean = np.mean(total_flux)

    FWHM = (max(fwhm_flux)+mean)/2.0    #The actual FWHM flux value
    
    print hdulist[i] 

    FWHM_y = []
    FWHM_x = []
    for j in range(len(total_flux)):
#        if total_flux[j] >= FWHM:
            if x[j] > l+7 and x[j] < NII_2-7:    #Only place values in FWHM_y and FWHM_x between the two NII lines (-10 from NII_2 as a quick fix for 91924 source)
                FWHM_y.append(total_flux[j])
                FWHM_x.append(x[j])

    hrange = [] #list of True/False values. True for wavelength slices inside of FWHM range. False otherwise. To be used to choose slices to include in collapsed cube.
    
    for j in range(len(x)):
        if x[j] < FWHM_x[0]:
            hrange.append(False)
        elif x[j] < FWHM_x[len(FWHM_x)-1]:
            hrange.append(True)
        else:
            hrange.append(False)

    summed_data = np.zeros([50,50])
    data = hdulist1[0].data #Doing this again to bring back the NaN values

    for j in range(0,len(data)):
        if hrange[j] == True:
            summed_data += data[j]  #summed_data is summed data within FWHM profile. ie. the collapsed Halpha image

    hdulist1[0].data = summed_data

    hdulist1.writeto(filename_list[i]+'collapsed.fits')

    #Plot
    line1, = plt.plot(x,total_flux)
    line5, = plt.plot((x[0],x[2047]),(0.0,0.0))
    line6, = plt.plot((x[0],x[2047]),(mean,mean))
    line7, = plt.plot((x[0],x[2047]),(FWHM,FWHM))
    line8, = plt.plot((FWHM_x[0],FWHM_x[0]),(min(total_flux),max(total_flux)))
    line9, = plt.plot((FWHM_x[len(FWHM_x)-1],FWHM_x[len(FWHM_x)-1]),(min(total_flux),max(total_flux)))
    line2, = plt.plot((l,l),(min(total_flux),max(total_flux)))
    line3, = plt.plot((NII_1,NII_1),(min(total_flux),max(total_flux)))
    line4, = plt.plot((NII_2,NII_2),(min(total_flux),max(total_flux)))
    line2.set_dashes([4,2,2,2])
    line3.set_dashes([2,2])
    line4.set_dashes([2,2])
    plt.title('Total Flux vs Wavelength from',fontsize=10)
    plt.xlabel('Wavelength (microns)',fontsize=10)
    plt.ylabel('Total Flux (ergs cm$\mathregular{^{-2}}$ $\\AA\mathregular{^{-1}}$ s$\mathregular{^{-1}}$)',fontsize=10)
    plt.axis('tight')
    plt.xlim([NII_1,NII_2])
    plt.ylim([0,max(fwhm_flux)])
    fig = plt.gcf()
    fig.set_size_inches(18.5, 5.5)
    plt.setp(line1, linewidth=0.1, color='b')
    plt.setp(line2, linewidth=0.1, color='g')
    plt.setp(line3, linewidth=0.1, color='r')
    plt.setp(line4, linewidth=0.1, color='r')
    plt.setp(line5, linewidth=0.1)
    plt.setp(line6, linewidth=0.1)
    plt.setp(line7, linewidth=0.1)
    plt.setp(line8, linewidth=0.1)
    plt.setp(line9, linewidth=0.1)
    plt.savefig(filename_list[i]+'collapsed.pdf')
    plt.close()

    hdulist1.close()
