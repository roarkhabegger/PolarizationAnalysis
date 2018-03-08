# -*- coding: utf-8 -*-
"""
Created on Sun Dec 24 12:27:25 2017

@author: Roark Habegger
"""
import glob
import numpy as np
import astropy.io.fits as apfits
import matplotlib.pyplot as plt
#from mpl_toolkits.mplot3d import Axes3D

#Open Server images
Serverfiles = glob.glob("C:\\Users\\Roark Habegger\\Documents\\PromptServer_.fits\\*.fits")
Serverimages = np.zeros((len(Serverfiles),2048,2048))

c=0
for i in Serverfiles:
	hdulist=apfits.open(i)
	#print hdulist[0].data.shape
	Serverimages[c]= np.array(hdulist[0].data)
	c+=1
print 'Server Images Open'

#open client images
Clientfiles = glob.glob("C:\\Users\\Roark Habegger\\Documents\\PromptClient_.fits\\*.fits")
Clientimages = np.zeros((len(Clientfiles),1024,1024))

c=0
for i in Clientfiles:
	hdulist=apfits.open(i)
	#print hdulist[0].data.shape
	Clientimages[c]= np.array(hdulist[0].data)
	c+=1
print 'Client images open'

starPix = np.zeros_like(Serverimages)
for k in range(0,len(Serverimages)):
	print '*****image '+str(k)
	print np.median(Serverimages[k])
	print np.mean(Serverimages[k])
	numStarPix = 0
	c=0
	for row in Serverimages[k]:
		upperEnd = np.where(row>1400)[0]
		numStarPix+=len(upperEnd)

		for pixel in upperEnd:
			d = 0

			try:
				right = Serverimages[k][c][pixel+4]
				d=1
				left = Serverimages[k][c][pixel-4]
				if (right>=1400 or left>=1400):
					starPix[k][c][pixel]=Serverimages[k][c][pixel]
			except:
				if (d==0 and Serverimages[k][c][pixel-4]>=1400):
					starPix[k][c][pixel]=Serverimages[k][c][pixel]
				elif (d==1 and Serverimages[k][c][pixel+4]>=1400):
					starPix[k][c][pixel]=Serverimages[k][c][pixel]

		c+=1
	print numStarPix
	impPixels = np.array((np.where(starPix[k][:]!=0)[0],np.where(starPix[k][:]!=0)[1]))


	#plt.figure(k+1)
	#plt.clf()
	#plt.imshow(starPix[k])
	#plt.show()









