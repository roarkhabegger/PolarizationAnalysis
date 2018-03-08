import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

exampletheta = np.array([np.pi/8.0, (3.0*np.pi)/8.0, np.pi/2.0, (7.0*np.pi)/8.0])
exampleintensity = np.array([103.266 + 3, 190.706 - 2, 266.364 + 1, 207.294])

def examplefunc(thet):
    return 100 * np.cos(2.0 * (thet - 2.0)) + 200

def cosineregression(thetadata, intensitydata):
    thetadata = np.asarray(thetadata) #converts input to arrays if not already arrays
    intensitydata = np.asarray(intensitydata)
    def cosine(theta, a, c, d): #theta IN RADIANS. a is amplitude, c is offset in intensity, d is shift in theta, b is 2pi/period 
        return a * np.cos(2.0 * (theta - d)) + c #np.cos IN RADIANS
    aguess = (np.amax(intensitydata) - np.amin(intensitydata)) / 2.0
    #bguess = 2.0 #(2*pi) / pi, same period for any data for our usage
    cguess = np.average(intensitydata)
    dguess = thetadata[np.argmax(intensitydata)] #argmax gives INDEX of max val
    paramsguess = [aguess, cguess, dguess]#initial guess for parameters
    optimized = curve_fit(cosine, exampletheta, exampleintensity, paramsguess)
    afit = optimized[0][0]
    cfit = optimized[0][1]
    dfit = optimized[0][2]
    def fittedcosine(the):
        return afit * np.cos(2.0 * (the - dfit)) + cfit
    return fittedcosine



plt.plot(exampletheta, exampleintensity, 'k.')


thetarange = np.arange(0, np.pi, 0.01)
plt.plot(thetarange, examplefunc(thetarange), 'r') #true function
fittedtest = cosineregression(exampletheta, exampleintensity)
plt.plot(thetarange, fittedtest(thetarange), 'b') #fitted function