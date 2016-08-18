# Fit exponential curve to data with convolution. Ydata will be actual radial profile values. Second array of convolution will be PSF. curve_fit from scipy.optimize is used to fit model to data.

import numpy as np
from scipy.optimize import curve_fit
from astropy.io import fits
import os
import fnmatch
import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
from astropy.convolution import convolve, convolve_fft, Gaussian1DKernel

gauss_kernel = Gaussian1DKernel(stddev=2)

# convolved model
def func(R, n, I_0, R_e):
    '''
    I_R = []
    for i in range(len(R)):
        I_R.append(I_0 * np.e**(-(1.9992*n-0.3721)*(R[i]/R_e)**(1.0/n)))
    return convolve_fft(I_R,np.array([5,4.3,2.1,1.0,0.3,0.2,0.1,0.01]))
    '''
    return convolve(I_0 * np.e**(-(1.9992*n-0.3721)*(R/R_e)**(1.0/n)), gauss_kernel)   #second array should be 1D PSF

# unconvolved model, maybe to use to print right parameters? Or maybe func does print right parameters...
def func2(R, n, I_0, R_e):
    return I_0 * np.e**(-(1.9992*n-0.3721)*(R/R_e)**(1.0/n))

Rdata = np.linspace(0,8,8)  #this will be R values of radial profile of target

#Y = func(Rdata, 1.0, 5, 12)
#Ydata = Y + 0.2 * np.random.normal(size=len(Rdata))

Ydata = np.array([5,4.3,2.1,1.0,0.3,0.2,0.1,0.01])  #this will be radial profile of target

curveFit = curve_fit(func, Rdata, Ydata, bounds = (0.001, [np.inf,np.inf,np.inf]))
print curveFit
curveFit2 = curve_fit(func2, Rdata, Ydata, bounds = (0.001, [np.inf,np.inf,np.inf]))
print curveFit2

Y = func(Rdata, curveFit[0][0], curveFit[0][1], curveFit[0][2])
Y2 = func2(Rdata, curveFit2[0][0], curveFit2[0][1], curveFit2[0][2])

plt.plot(Rdata, Y)
plt.plot(Rdata, Y2)
plt.show()
