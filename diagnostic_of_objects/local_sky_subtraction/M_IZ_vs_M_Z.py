#!/usr/bin/python

####
# Script that will plot the calculated M_IZ magnitudes from local_sky_subtraction.py vs M_Z magntiudes from the txt table files
####

from modules import *
from scipy import stats
from local_sky_subtraction import m_IZ_list, m_Z_list

#print m_IZ_list
#print m_Z_list

m_IZ_list = np.array(m_IZ_list)
m_Z_list = np.array(m_Z_list)

slope, intercept, r_value, p_value, std_err = stats.linregress(m_Z_list, m_IZ_list)
r_squared = r_value**2.0
best_fit = slope * m_Z_list + intercept

print r_squared

plt.plot(m_Z_list, m_IZ_list, 'o')
line1, = plt.plot(m_Z_list, best_fit)
plt.setp(line1, linewidth=0.1)
plt.title('Calculated IZ Magnitude vs Expected Z Magnitude')
plt.xlabel('Expected Z Magnitude')
plt.ylabel('Calculated IZ Magnitude')
plt.savefig('./figures/M_IZ_vs_M_Z.pdf')
plt.close()
