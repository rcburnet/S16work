#!/usr/bin/python

####
# Script that will look into combined OB data cubes found in /home/rburnet/reflex/project_data/OBs_combined/sky_subtracted/ , do a local sky subtraction, calculate magnitude of target, and plot the IZ and YJ spectrum of each target as detailed in the email sent by Dr. Balogh on June 15. Note: Plotted spectra are not total flux, but median flux for each wavelength slice. Magnitude is calculated by taking the sum of the aperture flux at each wavelength slice and choosing the median of all those wavelength slices as the flux and then calculating mag = -2.5 * log(flux/0.1*F0)
####

#Import modules
from modules import *
from spaxel_define import *

#Create list of filenames
hdulist=[]

for root,dirnames,filenames in os.walk('/home/rburnet/reflex/project_data/OBs_combined/sky_subtracted/'):
        for filename in fnmatch.filter(filenames, 'COMBINE_SCI_RECONSTRUCTED_[0-9]*'):
                hdulist.append(os.path.join(root, filename))

hdulist = sorted(hdulist)

spaxel_array_select = np.array([CL0034_IZ_Target_10, CL0034_IZ_Target_11, CL0034_IZ_Target_12, CL0034_IZ_Target_14, CL0034_IZ_Target_16, CL0034_IZ_Target_2, CL0034_IZ_Target_20, CL0034_IZ_Target_21, CL0034_IZ_Target_22, CL0034_IZ_Target_23, CL0034_IZ_Target_24, CL0034_IZ_Target_25, CL0034_IZ_Target_26, CL0034_IZ_Target_30, CL0034_IZ_Target_31, CL0034_IZ_Target_4, CL0034_IZ_Target_5, CL0034_IZ_Target_6, CL0034_IZ_Target_7, CL0034_IZ_Target_8, CL0034_YJ_Target_10, CL0034_YJ_Target_11, CL0034_YJ_Target_12, CL0034_YJ_Target_14, CL0034_YJ_Target_16, CL0034_YJ_Target_2, CL0034_YJ_Target_20, CL0034_YJ_Target_21, CL0034_YJ_Target_22, CL0034_YJ_Target_23, CL0034_YJ_Target_24, CL0034_YJ_Target_25, CL0034_YJ_Target_26, CL0034_YJ_Target_31, CL0034_YJ_Target_4, CL0034_YJ_Target_5, CL0034_YJ_Target_6, CL0034_YJ_Target_7, CL0034_YJ_Target_8, CL0036_IZ_Target_10, CL0036_IZ_Target_11, CL0036_IZ_Target_12, CL0036_IZ_Target_13, CL0036_IZ_Target_15, CL0036_IZ_Target_17, CL0036_IZ_Target_19, CL0036_IZ_Target_22, CL0036_IZ_Target_26, CL0036_IZ_Target_28, CL0036_IZ_Target_3, CL0036_IZ_Target_30, CL0036_IZ_Target_4, CL0036_IZ_Target_5, CL0036_IZ_Target_6, CL0036_IZ_Target_7, CL0036_IZ_Target_8, CL0036_YJ_Target_10, CL0036_YJ_Target_11, CL0036_YJ_Target_12, CL0036_YJ_Target_13, CL0036_YJ_Target_15, CL0036_YJ_Target_17, CL0036_YJ_Target_18, CL0036_YJ_Target_19, CL0036_YJ_Target_28, CL0036_YJ_Target_3, CL0036_YJ_Target_30, CL0036_YJ_Target_34, CL0036_YJ_Target_4, CL0036_YJ_Target_5, CL0036_YJ_Target_6, CL0036_YJ_Target_7, CL0036_YJ_Target_8])   #all the arrays from spaxel_define that define the apertures for each target

### Read txt table files for target/arm info
file_CL0034 = open('/home/rburnet/S16work/txt_tables/KMOS_GCLASS_CLUS0034.txt','r')
file_CL0036_YJ = open('/home/rburnet/S16work/txt_tables/KMOS_GCLASS_CLUS0036_YJ.txt','r')
file_CL0036_IZ = open('/home/rburnet/S16work/txt_tables/KMOS_GCLASS_CLUS0036_IZ.txt','r')

CL0034_info = file_CL0034.readlines()
CL0036_YJ_info = file_CL0036_YJ.readlines()
CL0036_IZ_info = file_CL0036_IZ.readlines()

file_CL0034.close()
file_CL0036_YJ.close()
file_CL0036_IZ.close()

for i in range(len(CL0034_info)):
    CL0034_info[i] = CL0034_info[i].split(' ')
    for j in range(len(CL0034_info[i])):
        if j < 3:
            CL0034_info[i][j] = int(CL0034_info[i][j])
        else:
            CL0034_info[i][j] = float(CL0034_info[i][j])

for i in range(len(CL0036_YJ_info)):
    CL0036_YJ_info[i] = CL0036_YJ_info[i].split(' ')
    for j in range(len(CL0036_YJ_info[i])):
        if j < 3:
            CL0036_YJ_info[i][j] = int(CL0036_YJ_info[i][j])
        else:
            CL0036_YJ_info[i][j] = float(CL0036_YJ_info[i][j])

for i in range(len(CL0036_IZ_info)):
    CL0036_IZ_info[i] = CL0036_IZ_info[i].split(' ')
    for j in range(len(CL0036_IZ_info[i])):
        if j < 3:
            CL0036_IZ_info[i][j] = int(CL0036_IZ_info[i][j])
        else:
            CL0036_IZ_info[i][j] = float(CL0036_IZ_info[i][j])

#Read data cubes and do calculations/plotting

x_YJ = np.linspace(1.0,1.359,2048)
x_IZ = np.linspace(0.78,1.09,2048)

for i in range(len(hdulist)):
    if not os.path.isfile('./figures/'+hdulist[i][62:]+'.pdf'):
        hdulist1 = fits.open(hdulist[i])
        target_flux = []
        median_flux = []
        sky_flux = []        

        #Loop into each wavelength slice of every data cube and create specified arrays from data
        for j in range(len(hdulist1[1].data)):
            spaxel_array = hdulist1[1].data[j]
            spaxel_array_target = spaxel_array[spaxel_array_select[i]]
            spaxel_array_sky = spaxel_array[~spaxel_array_select[i]]
            spaxel_array_sky = spaxel_array_sky[~np.isnan(spaxel_array_sky)]
            median_sky_flux = np.median(spaxel_array_sky)
            spaxel_array = spaxel_array - median_sky_flux
            spaxel_array_target = spaxel_array_target - median_sky_flux
            spaxel_array_target = spaxel_array_target[~np.isnan(spaxel_array_target)]
            spaxel_array_no_NaN = spaxel_array[~np.isnan(spaxel_array)] #originally had a total_flux list that would be the sum of this array, not sure if I still need it, kept it here anyways just in case I'll need it sometime in the future
            target_flux.append(np.sum(spaxel_array_target)) #to be used in magnitude calculation
            median_flux.append(np.median(spaxel_array_target))  #to be used in plotting spectrum
            sky_flux.append(median_sky_flux)    #to be used in plotting spectrum for local sky
        
        target_flux = [x for x in target_flux if x != target_flux[0]]   #ignore values of 0.0 flux. eg. beginning and end of data
        total_target_flux = np.median(target_flux)
        
        #Magnitude calculation. F0 is zero point flux for the filter. total_target_flux is the median of all wavelength slice's aperture's summed flux.
        if 'YJ' in hdulist[i]:
            F0 = 3.129e-9
        if 'IZ' in hdulist[i]:
            F0 = 7.63e-9
        mag = -2.5*np.log10(total_target_flux/(0.1*F0))
        #print hdulist[i][62:], mag

        #Now plot spectrum

        if not os.path.exists('./figures/'+hdulist[i][62:72]):
            os.makedirs('./figures/'+hdulist[i][62:72])

        if 'CL0034' in hdulist[i]:
            for k in range(len(CL0034_info)):
                if '_' in hdulist[i][-7:-5]:
                    if CL0034_info[k][0] == int(hdulist[i][-6:-5]):
                        z = CL0034_info[k][5]
                else:
                    if CL0034_info[k][0] == int(hdulist[i][-7:-5]):
                        z = CL0034_info[k][5]
        if 'CL0036' in hdulist[i]:
            if 'IZ' in hdulist[i]:
                for k in range(len(CL0036_IZ_info)):
                    if '_' in hdulist[i][-7:-5]:
                        if CL0036_IZ_info[k][0] == int(hdulist[i][-6:-5]):
                            z = CL0036_IZ_info[k][5]
                    else:
                        if CL0036_IZ_info[k][0] == int(hdulist[i][-7:-5]):
                            z = CL0036_IZ_info[k][5]
        if 'CL0036' in hdulist[i]:
            if 'YJ' in hdulist[i]:
                for k in range(len(CL0036_YJ_info)):
                    if '_' in hdulist[i][-7:-5]:
                        if CL0036_YJ_info[k][0] == int(hdulist[i][-6:-5]):
                            z = CL0036_YJ_info[k][5]
                    else:
                        if CL0036_YJ_info[k][0] == int(hdulist[i][-7:-5]):
                            z = CL0036_YJ_info[k][5]
        #print hdulist[i], z
        if 'YJ' in hdulist[i]:
            #Plot Halpha
            line1, = plt.plot(x_YJ,median_flux)
            line6, = plt.plot(x_YJ,sky_flux)
            minimum = min([min(median_flux),min(sky_flux)])
            maximum = max([max(median_flux),max(sky_flux)])
            y_range = (minimum,maximum)
            l = 0.6563 * ( 1 + z )
            NII_1 = 0.6548 * ( 1 + z )
            NII_2 = 0.6584 * ( 1 + z )
            line2, = plt.plot((l,l),y_range)
            line3, = plt.plot((NII_1,NII_1),(-2e-19,3e-19))
            line4, = plt.plot((NII_2,NII_2),(-2e-19,3e-19))
            line5, = plt.plot((1.0,1.35882),(0.0,0.0))
            line2.set_dashes([4,2,2,2])
            line3.set_dashes([2,2])
            line4.set_dashes([2,2])
            plt.title('Total Flux vs Wavelength centered on Halpha from \n'+hdulist[i],fontsize=10)
            plt.xlabel('Wavelength (microns)',fontsize=10)
            plt.ylabel('Total Flux (ergs cm$\mathregular{^{-2}}$ $\\AA\mathregular{^{-1}}$ s$\mathregular{^{-1}}$)',fontsize=10)
            plt.axis('tight')

            # pick minimum and maximum y values
            flux_domain = []
            sky_domain = []
            for m in range(len(x_IZ)):
                if x_YJ[m] > l-0.02 and x_YJ[m] < l+0.02:
                    flux_domain.append(median_flux)
                    sky_domain.append(sky_flux)
            minimum = min(np.mean(flux_domain),np.mean(sky_domain))
            maximum = max(np.mean(flux_domain),np.mean(sky_domain))
            plt.ylim([minimum-2e-19,maximum+2e-19])

            plt.xlim([l-0.02,l+0.02])
            fig = plt.gcf()
            fig.set_size_inches(18.5, 5.5)
            plt.setp(line1, linewidth=0.1, color='b')
            plt.setp(line2, linewidth=0.1, color='g')
            plt.setp(line5, linewidth=0.1, color='k')
            plt.setp(line6, linewidth=0.1, color='y')
            print './figures/'+hdulist[i][62:-5]+'_Halpha.pdf'
            plt.savefig('./figures/'+hdulist[i][62:-5]+'_Halpha.pdf')
            plt.close()
            hdulist1.close()

        if 'IZ' in hdulist[i]:
            #Plot Hbeta
            line1, = plt.plot(x_IZ,median_flux)
            line6, = plt.plot(x_IZ,sky_flux)
            minimum = min([min(median_flux),min(sky_flux)])
            maximum = max([max(median_flux),max(sky_flux)])
            y_range = (minimum,maximum)
            Hbeta = 0.4861 * ( 1 + z )
            line2, = plt.plot((Hbeta,Hbeta), (-6e-19,6e-19))
            line2.set_dashes([4,2,2,2])
            line5, = plt.plot((0.78,1.08985),(0.0,0.0))
            plt.title('Total Flux vs Wavelength centered on Hbeta from \n'+hdulist[i],fontsize=10)
            plt.xlabel('Wavelength (microns)',fontsize=10)
            plt.ylabel('Flux (ergs cm$\mathregular{^{-2}}$ $\\AA\mathregular{^{-1}}$ s$\mathregular{^{-1}}$)',fontsize=10)
            plt.axis('tight')

            flux_domain = []
            sky_domain = []
            for m in range(len(x_IZ)):
                if x_IZ[m] > Hbeta-0.02 and x_IZ[m] < Hbeta+0.02:
                    flux_domain.append(median_flux)
                    sky_domain.append(sky_flux)
            minimum = min(np.mean(flux_domain),np.mean(sky_domain))
            maximum = max(np.mean(flux_domain),np.mean(sky_domain))
            plt.ylim([minimum-2e-19,maximum+2e-19])

            plt.xlim([Hbeta-0.02,Hbeta+0.02])
            plt.setp(line1, linewidth=0.1, color='b')
            plt.setp(line2, linewidth=0.1, color='g' )
            plt.setp(line5, linewidth=0.1, color='k')
            plt.setp(line6, linewidth=0.1, color='y')
            fig = plt.gcf()
            fig.set_size_inches(18.5, 5.5)
            print './figures/'+hdulist[i][62:-5]+'_Hbeta.pdf'
            plt.savefig('./figures/'+hdulist[i][62:-5]+'_Hbeta.pdf')
            plt.close()

            #Plot OIII
            line1, = plt.plot(x_IZ,median_flux)
            line6, = plt.plot(x_IZ,sky_flux)
            minimum = min([min(median_flux),min(sky_flux)])
            maximum = max([max(median_flux),max(sky_flux)])
            y_range = (minimum,maximum)
            OIII = 0.5007 * ( 1 + z )
            line3, = plt.plot((OIII,OIII), (-6e-19,6e-19))
            line3.set_dashes([4,2,2,2])
            line5, = plt.plot((0.78,1.08985),(0.0,0.0))
            plt.title('Total Flux vs Wavelength centered on Hbeta from \n'+hdulist[i],fontsize=10)
            plt.xlabel('Wavelength (microns)',fontsize=10)
            plt.ylabel('Flux (ergs cm$\mathregular{^{-2}}$ $\\AA\mathregular{^{-1}}$ s$\mathregular{^{-1}}$)',fontsize=10)
            plt.axis('tight')

            flux_domain = []
            sky_domain = []
            for m in range(len(x_IZ)):
                if x_IZ[m] > OIII-0.02 and x_IZ[m] < OIII+0.02:
                    flux_domain.append(median_flux)
                    sky_domain.append(sky_flux)
            minimum = min(np.mean(flux_domain),np.mean(sky_domain))
            maximum = max(np.mean(flux_domain),np.mean(sky_domain))
            plt.ylim([minimum-2e-19,maximum+2e-19])

            plt.xlim([OIII-0.02,OIII+0.02])
            plt.setp(line1, linewidth=0.1, color='b')
            plt.setp(line3, linewidth=0.1, color='r')
            plt.setp(line5, linewidth=0.1, color='k')
            plt.setp(line6, linewidth=0.1, color='y')
            fig = plt.gcf()
            fig.set_size_inches(18.5, 5.5)
            print './figures/'+hdulist[i][62:-5]+'_OIII.pdf'
            plt.savefig('./figures/'+hdulist[i][62:-5]+'_OIII.pdf')
            plt.close()
            hdulist1.close()
