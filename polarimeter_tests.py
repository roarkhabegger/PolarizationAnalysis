import polarimeteranalysis_p36 as pol

import numpy as np
import matplotlib.pyplot as plt

a = 1
b = 0.5
d = 1

def ellipse(theta):
    return (a*b) / np.sqrt((b*np.cos(theta-d))**2.0 + (a*np.sin(theta-d))**2.0)

def ellipsePlot(fluxdata): # plots an ellipse from flux data
    angles = np.arange(0, 2.0*np.pi + np.pi/8.0, np.pi/8.0)
    fluxes = []
    for j in range(2):
        for i in range(0,4,2):
            fluxes.append(fluxdata[i][j])
        for i in range(1,5,2):
            fluxes.append(fluxdata[i][j])

    #print(fluxes)

    fluxes += fluxes
    fluxes.append(fluxes[0])

    plt.subplot(111, projection='polar')
    plt.plot(angles, fluxes)
    plt.grid(True)
    plt.title("Flux vs. angle")
    plt.show()


def testAnalysis(fluxdata, sigmadata):
    results = pol.polarimetry(fluxdata, sigmadata)


    print("\n\nPolarization 1 = " + str(results[0][0]) + " +/- " + str(results[0][2]) + "\n" 
          "Polarization Angle 1 = " + str(results[0][1]) + " +/- " + str(results[0][3]) + "\n" 
          "Polarization 2 = " + str(results[1][0]) + " +/- " + str(results[1][2]) + "\n" 
          "Polarization Angle 2 = " + str(results[1][1]) + " +/- " + str(results[0][3]) + "\n" 
          "Average Polarization = " + str(results[2][0]) + " +/- " + str(results[2][2]) + "\n" 
          "Average Polarization Angle = " + str(results[2][1]) + " +/- " + str(results[2][3])


          )
    ellipsePlot(fluxdata)

# test sigma data
sigmadata = ((0.05, 0.1), (0.08, 0.07), (0.06, 0.09), (0.1, 0.04))

# Circular Case: no polarization
fluxdata_circ = ((2.0, 2.0) , (2.0,2.0) , (2.0,2.0) , (2.0,2.0))

# Linear Case: complete polarization
fluxdata_linear = ((2.0, 0.0) , (0.0, 0.0), (0.0, 0.0), (0.0, 0.0))

# Linear Case shifted up with polarization angle of 22.5
fluxdata_linear_shifted = ((0.0, 0.0) , (0.0, 0.0), (2.0, 0.0), (0.0, 0.0))

#Semi polarized along axis
fluxdata_semi = ((1.0,0.5) , (0.634256,0.634256), (0.833524,0.52995), (0.52995, 0.833524))

#Same, but shifted up by d radians
d = -1
fluxdata_semi_shifted = [[0.0,0.0], [0.0,0.0], [0.0,0.0], [0.0,0.0]]

c = 0
for j in range(2):
        for i in range(0,4,2):
            fluxdata_semi_shifted[i][j] = ellipse(c)
            c += np.pi/8.0
        for i in range(1,5,2):
            fluxdata_semi_shifted[i][j] = ellipse(c)
            c += np.pi/8.0
# main

testAnalysis(fluxdata_semi_shifted, sigmadata)