#import imexam
import numpy as np
import matplotlib.pyplot as plt
from astropy.modeling.models import Sersic1D

x0 = np.linspace(0,15,87)

def S(R):
    A = 1.428430e-01
    B = -4.265817e-02
    C = 5.192784e-03
    D = -2.909510e-04
    E = 6.164860e-06
    F = -5.493876e-02
    Sigma = 6.700000e-01
    return A + B * R + C * R**2.0 + D * R**3.0 + E * R**4.0 + F * np.e**(-R/(2 * Sigma**2.0))

y = S(x0)

cumul = np.cumsum(y)/max(np.cumsum(y))    #cumulative sum of radial profile
'''
imexam.list_active_ds9()    #list active ds9 sessions. Must have ds9 open.
viewer=imexam.connect('/tmp/.xpa/DS9_ds9.5228') #This should be the output of the above line. Need to change it to that.
viewer.load_fits('/home/rburnet/SAMI/data/without_HI_detections/77710_red_7_Y13SAR1_P014_15T029.fits')
viewer.imexam()
'''
x = np.array([0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15])
y = np.array([0.000000, 0.440216, 1.090727, 2.089978, 2.997073, 4.029789, 4.736484, 5.309045, 5.938111, 6.469439, 6.921639, 7.237229, 7.503949, 7.793467, 8.004063, 8.085778])

y = y/max(y)

radialProfile = open('/home/rburnet/.AperturePhotometryTool/radialProfile.dat')
read = radialProfile.readlines()
x1 = []
y1 = []
for i in range(len(read)):
    read[i] = read[i].split(' ')
    x1.append(float(read[i][0]))
    y1.append(float(read[i][5]))

y1 = np.array(y1)
cumul1 = np.cumsum(y1)/max(np.cumsum(y1))

x2 = np.linspace(0,15,10000)
y2 = S(x2)
cumul2 = np.cumsum(y2)/max(np.cumsum(y2))

#plt.plot(x0,cumul)
#plt.plot(x2,cumul2)
plt.plot(x,y)
plt.plot(x1,cumul1)
#plt.plot(x1,cumul1)
plt.show()

