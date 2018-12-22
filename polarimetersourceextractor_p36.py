# -*- coding: utf-8 -*-
"""
Created on Sat May 19 18:29:43 2018

@author: Nick Konz
"""
import sep
import numpy as np
import astropy.io.fits as fits
import matplotlib.pyplot as plt
from matplotlib.patches import Ellipse

# read image into standard 2d np array:
filename1 = 'Server2.fits'
filename2 = 'client2.fits'

hdul1 = fits.open(filename1) #list-like collection of HDU objects
data1 = hdul1[1].data
data1 = data1.byteswap().newbyteorder() #byteswapping is necessary so that background subtraction works

hdul2 = fits.open(filename2) #list-like collection of HDU objects
data2 = hdul2[1].data
data2 = data2.byteswap().newbyteorder() #byteswapping is necessary so that background subtraction works

#TEST PLOTTING
'''
plt.imshow(data, cmap= 'gray')
plt.colorbar();
'''

#data must be background subtracted before sources can be detected
#in SEP, background estimation and source detection are two seperate steps

bkg1 = sep.Background(data1, mask=None, bw=64, bh=64, fw=3, fh=3)
bkg2 = sep.Background(data2, mask=None, bw=64, bh=64, fw=3, fh=3)
'''
# get a "global" mean and noise of the image background:
print(bkg.globalback)
print(bkg.globalrms)
'''

'''
# evaluate background as 2-d array, same size as original image
bkg_image = bkg.back()
# bkg_image = np.array(bkg) # equivalent to above

# show the background
plt.imshow(bkg_image, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
'''

'''
# evaluate the background noise as 2-d array, same size as original image
bkg_rms = bkg.rms()

# show the background noise
plt.imshow(bkg_rms, interpolation='nearest', cmap='gray', origin='lower')
plt.colorbar();
'''

# subtract the background
data_sub1 = data1 - bkg1
data_sub2 = data2 - bkg2

#now that background is subtracted, can run object detection
#here we use the detection threshold to be threshold times the global background RMS
thresh = 70
objects1 = sep.extract(data_sub1, thresh, err=bkg1.globalrms)
objects2 = sep.extract(data_sub2, thresh, err=bkg2.globalrms)

#how many objects were detected
print(len(objects1))
print(len(objects2))

#aperture photometry test
#example of simple circular aperture photometry 
#with a 3 pixel radius at the locations of the objects:

'''Objects 1'''
flux1, fluxerr1, flag1 = sep.sum_circle(data_sub1, objects1['x'], objects1['y'],
                                     3.0, err=bkg1.globalrms, gain=1.0)

# show the results:
for i in range(len(objects1)):
    print("object {:d}: flux1 = {:f} +/- {:f}".format(i, flux1[i], fluxerr1[i]))
    
# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(data_sub1), np.std(data_sub1)
im = ax.imshow(data_sub1, cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(objects1)):
    e = Ellipse(xy=(objects1['x'][i], objects1['y'][i]),
                width=6*objects1['a'][i],
                height=6*objects1['b'][i],
                angle=objects1['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)

'''Objects 2'''
flux2, fluxerr2, flag2 = sep.sum_circle(data_sub2, objects2['x'], objects2['y'],
                                     3.0, err=bkg2.globalrms, gain=1.0)

# show the results:
for i in range(len(objects2)):
    print("object {:d}: flux2 = {:f} +/- {:f}".format(i, flux2[i], fluxerr2[i]))
    
# plot background-subtracted image
fig, ax = plt.subplots()
m, s = np.mean(data_sub2), np.std(data_sub2)
im = ax.imshow(data_sub2, cmap='gray',
               vmin=m-s, vmax=m+s, origin='lower')

# plot an ellipse for each object
for i in range(len(objects2)):
    e = Ellipse(xy=(objects2['x'][i], objects2['y'][i]),
                width=6*objects2['a'][i],
                height=6*objects2['b'][i],
                angle=objects2['theta'][i] * 180. / np.pi)
    e.set_facecolor('none')
    e.set_edgecolor('red')
    ax.add_artist(e)
