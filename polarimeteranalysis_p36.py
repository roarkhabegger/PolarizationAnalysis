'''
@author: Nick Konz
'''
import numpy as np

def polarizationUncert(p,Theta,dataset,sigmaset): #assuming dataset is ((F0,F90),(F45,F135)) (or plus 22.5). sigmas organized same way
    F0 = dataset[0][0]
    F90 = dataset[0][1]
    F45 = dataset[1][0]
    F135 = dataset[1][1]

    sig0 = sigmaset[0][0]
    sig90 = sigmaset[0][1]
    sig45 = sigmaset[1][0]
    sig135 = sigmaset[1][1]

    sigma_p = np.sqrt(#adding up the various partials in quadrature (standard propagation of error)
        sig0**2.0 * ( -1.0*np.sqrt((F45-F135)**2.0 + (F0-F90)**2.0) / (F0+F90)**2.0 + (F0-F90) / ( np.sqrt((F45-F135)**2.0 + (F0-F90)**2.0) * (F0+F90) ) )**2.0 +
        sig90**2.0 * ( -1.0*np.sqrt((F45-F135)**2.0 + (F0-F90)**2.0) / (F0+F90)**2.0 - (F0-F90) / ( np.sqrt((F45-F135)**2.0 + (F0-F90)**2.0) * (F0+F90) ) )**2.0 + #this aren't EXACTLY the same... see the minus sign
        sig45**2.0 * ((F45-F135) / ( np.sqrt((F45-F135)**2.0 + (F0-F90)**2.0) * (F0+F90) ))**2.0 +
        sig135**2.0 * (-1.0*(F45-F135) / ( np.sqrt((F45-F135)**2.0 + (F0-F90)**2.0) * (F0+F90) ))**2.0
        )

    sigma_Theta = np.sqrt(
        sig0**2.0 * ((F135-F45) / (2.0*(1.0 + (F45-F135)**2.0/(F0-F90)**2.0)*(F0-F90)**2.0))**2.0 +
        sig90**2.0 * ((F45-F135) / (2.0*(1.0 + (F45-F135)**2.0/(F0-F90)**2.0)*(F0-F90)**2.0))**2.0 +
        sig45**2.0 * (1.0 / (2.0*(1.0 + (F45-F135)**2.0/(F0-F90)**2.0)*(F0-F90)))**2.0 +
        sig135**2.0 * (-1.0 / (2.0*(1.0 + (F45-F135)**2.0/(F0-F90)**2.0)*(F0-F90)))**2.0
        )


    return sigma_p,sigma_Theta


def polarizationParams(dataset):

    F0 = dataset[0][0]
    F90 = dataset[0][1]
    F45 = dataset[1][0]
    F135 = dataset[1][1]

    p = np.sqrt( (F0-F90)**2.0 + (F45-F135)**2.0 ) / (F0 + F90)
    Theta = 0.5 * np.arctan( (F45 - F135) / (F0-F90) )

    return p, Theta

def polarimetry(data,sigmadata): #assumes that data is a list of ordered fluxes by angle for each camera pair, i.e. for angles of ((0,90),(22.5,112.5),(45,135),(67.5,157.5))

    #returns a tuple of ((p, sigma_p), (Theta, sigma_Theta))


    dataset1 = (data[0], data[2]) #((F0,F90) , (F45,F135))
    dataset2 = (data[1], data[3]) #((F225,F1125) , (F675, F1575)) #not sure about the validity of this calculation

    sigmaset1 = (sigmadata[0], sigmadata[2])
    sigmaset2 = (sigmadata[1], sigmadata[3])

    p1, Theta1 = polarizationParams(dataset1)
    p2, Theta2 = polarizationParams(dataset2)

    p_avg = (p1 + p2) / 2.0
    Theta_avg = (Theta1 + Theta2) / 2.0

    #Uncertainty Calculation:
    #Using standard propagation of uncertainty:

    sigma_p1, sigma_Theta1 = polUncert(p1, Theta1, dataset1, sigmaset1)
    sigma_p2, sigma_Theta2 = polUncert(p2, Theta2, dataset2, sigmaset2)

    sigma_p_avg = 0.5 * np.sqrt(sigma_p1**2.0 + sigma_p2**2.0)
    sigma_Theta_avg = 0.5 * np.sqrt(sigma_Theta1**2.0 + sigma_Theta2**2.0)

    return ( (p1, Theta1, sigma_p1, sigma_Theta1), (p2, Theta2, sigma_p2, sigma_Theta2), (p_avg, Theta_avg, sigma_p_avg, sigma_Theta_avg) )