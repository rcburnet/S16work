#!/usr/bin/python

####
# Script that will plot the calculated M_IZ magnitudes from local_sky_subtraction.py vs M_Z magntiudes from the txt table files
####

from modules import *
from scipy import stats
from local_sky_subtraction import m_IZ_list, m_Z_list

#print m_IZ_list
#print m_Z_list

hduname = []

for i in range(len(m_IZ_list)):
    hduname.append(m_IZ_list[i][0])
    m_IZ_list[i] = m_IZ_list[i][1]

for i in range(len(hduname)):
    if '_' in hduname[i][-7:-5]:
        hduname[i] = hduname[i][-42:-36]+'_'+hduname[i][-6:-5]
    else:
        hduname[i] = hduname[i][-43:-37]+'_'+hduname[i][-7:-5]

m_IZ_list = np.array(m_IZ_list)
m_Z_list = np.array(m_Z_list)

slope, intercept, r_value, p_value, std_err = stats.linregress(m_Z_list, m_IZ_list)
r_squared = r_value**2.0
best_fit = slope * m_Z_list + intercept

print r_squared
fig = plt.figure()
ax = fig.add_subplot(111)
plt.scatter(m_Z_list, m_IZ_list, s=10)
k = 0
for i,j in zip(m_Z_list,m_IZ_list):
    if hduname[k] == 'CL0036_5' or hduname[k] == 'CL0036_22' or hduname[k] == 'CL0036_13' or hduname[k] == 'CL0036_7':
        ax.annotate('%s' % hduname[k], xy=(i,j), xytext=(0,2.5), textcoords='offset points', fontsize=3)
    else:
        ax.annotate('%s' % hduname[k], xy=(i,j), xytext=(2.5,0), textcoords='offset points', fontsize=3)
    k += 1
line1, = plt.plot(m_Z_list, best_fit)
ax.annotate('Slope = '+str('%.3f' % slope), xy=(max(m_Z_list), max(best_fit)), xytext=(0,2.5), textcoords='offset points', fontsize=5)
ax.annotate('r$^2$ = '+str('%.3f' % r_squared), xy=(max(m_Z_list), max(best_fit)), xytext=(0,-2.5), textcoords='offset points', fontsize=5)
plt.setp(line1, linewidth=0.1)
plt.title('Calculated IZ Magnitude vs Expected Z Magnitude')
plt.xlabel('Expected Z Magnitude')
plt.ylabel('Calculated IZ Magnitude')
plt.grid()
plt.savefig('./figures/M_IZ_vs_M_Z.pdf')
plt.close()
