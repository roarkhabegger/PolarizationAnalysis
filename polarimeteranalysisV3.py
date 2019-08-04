'''
@author: Nick Konz
'''
import numpy as np
from scipy.optimize import curve_fit
import matplotlib.pyplot as plt

ANGLES = ((0, np.pi/2.0),
          (np.pi/4.0, np.pi*0.75),
          (np.pi/8.0, np.pi*0.625),
          (np.pi*0.375, np.pi*0.875)
          )

def polRatio(theta, Pi, Theta):
    return Pi * np.cos(2.0*(theta - Theta))

def fitPol(angleData, ratioData, sigmaData):
    plt.figure(1)
    plt.plot(angleData, ratioData, "k*", markersize=12)

    plotX = np.linspace(np.amin(angleData), np.amax(angleData), 1000)
    plotFit = np.zeros(1000)


    fit = curve_fit(polRatio, angleData, ratioData, sigma=sigmaData, bounds=(np.array([0, 0]), np.array([np.inf, np.pi/2.0])))

    for i in range(plotFit.size):
        plotFit[i] = polRatio(plotX[i], fit[0][0], fit[0][1])

    plt.plot(plotX, plotFit, "-b")
    #plt.show()

    #note: uncertainties can only be determined for >1 degrees of freedom. Else, will return an inf.
    return(fit)

def polUncertainty(data, sigmadata):
    sigmaRatios = np.zeros(4)

    for i in range(4):
        F0 = data[i][0]
        F90 = data[i][1]

        s0 = sigmadata[i][0]
        s90 = sigmadata[i][1]

        sigmaRatios[i] = 2.0 * np.sqrt((s0 * (F90/((F0 + F90)**2.0)))**2.0 + (s90 * (F0/((F0 + F90)**2.0)))**2.0)

    return sigmaRatios


def polarimetry(data, sigmadata):
    #assumes that data is a list of ordered fluxes by angle for each camera pair, 
    #i.e. for angles of ((0,90),(45,135),(22.5,112.5),(67.5,157.5)). same for sigma data.

    sigmaRatios = polUncertainty(data, sigmadata) # size 4

    ratios = np.zeros(4)

    for j, ratio in enumerate(ratios):
        F0 = data[j][0]
        F90 = data[j][1]
        ratios[j] = (F0 - F90) / (F0 + F90)

    #returns a tuple of ((p, sigma_p), (Theta, sigma_Theta))

    allParams = [[],[]]
    allSigmas = [[],[]]

    #pairwise calculations                  (0, 45, 22.5, 67.5)
    for i in range(3):
        sigmaData = np.array([sigmaRatios[i], sigmaRatios[i+1]])
        ratioData = np.array([ratios[i], ratios[i+1]])
        angleData = np.array([ANGLES[i][0], ANGLES[i+1][0]])

        params, paramsCov = fitPol(angleData, ratioData, sigmaData)
    
        allParams[0].append(params[0])
        allParams[1].append(params[1])

        #uncertainty calculation
       
        allSigmas[0].append(paramsCov[0][0])
        allSigmas[1].append(paramsCov[1][1])

    #triplet calculations
    for i in range(2):
        sigmaData = np.array([sigmaRatios[i], sigmaRatios[i+1], sigmaRatios[i+2]])
        ratioData = np.array([ratios[i], ratios[i+1], ratios[i+2]])
        angleData = np.array([ANGLES[i][0], ANGLES[i+1][0], ANGLES[i+2][0]])

        params, paramsCov = fitPol(angleData, ratioData, sigmaData)
    
        allParams[0].append(params[0])
        allParams[1].append(params[1])

        #uncertainty calculation
       
        allSigmas[0].append(paramsCov[0][0])
        allSigmas[1].append(paramsCov[1][1])

    #quadruplet calculation

    #triplet calculations
    sigmaData = sigmaRatios
    ratioData = ratios
    angleData = np.array([ANGLES[0][0], ANGLES[1][0], ANGLES[2][0], ANGLES[3][0]])

    params, paramsCov = fitPol(angleData, ratioData, sigmaData)
    
    allParams[0].append(params[0])
    allParams[1].append(params[1])

    #uncertainty calculation
       
    allSigmas[0].append(paramsCov[0][0])
    allSigmas[1].append(paramsCov[1][1])


    #take average of everything

    p = allParams[0]
    Theta = allParams[1]

    sigma_p = allSigmas[0]
    sigma_Theta = allSigmas[1]
 
    return ((p, sigma_p), (Theta, sigma_Theta))