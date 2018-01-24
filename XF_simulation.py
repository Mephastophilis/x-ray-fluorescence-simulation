#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Thu Jan  4 14:01:35 2018

@author: bryanquigley
"""
from __future__ import division

import sys
sys.path.insert(0,"/Users/bryanquigley/Documents/Xray Luminescence and Fluorescence/XFXL Simulation")

import numpy as np
import matplotlib.pyplot as plt

#create capillary tube mask
tubemask = np.zeros((5000,2352), dtype=np.float64)
voxelsize =  0.0100 # 0.01 mm
#voxelsizex = 0.0459751172 #mm/pixel
voxelsizex = 0.0200 #mm
tubediam = 0.59436 # 0.6 mm
#tubediam = 0.8636 # 0.6 mm
spacing = 3.4 # 4 mm
tube1center = 1177 - (spacing/2)/voxelsizex - (tubediam/2)/voxelsizex
tube2center = 1177 + (spacing/2)/voxelsizex + (tubediam/2)/voxelsizex
Rsquared = (tubediam/2)

for i in xrange(5000):
    for j in xrange(2352):
        if j < 1176:

            if ((i-2500)*voxelsize) ** 2 + ((j-tube1center)*voxelsizex) ** 2 <= Rsquared:
                tubemask[i,j] = 1

        else:
            if ((i-2500)*voxelsize) ** 2 + ((j-tube2center)*voxelsizex) ** 2 <= Rsquared:
                tubemask[i,j] = 1

plt.imsave('tubemask.jpg', tubemask)
plt.imshow(tubemask)

#create gel phantom around tubes
phantommask = np.zeros((5000,2352), dtype=np.float64)
for j in xrange(5000):
    if j>=1700:
        phantommask[j,:]=1

for i in xrange(5000):
    for j in xrange(2352):
        if tubemask[i,j]==1:
            phantommask[i,j]=2

plt.imsave('phantommask.jpg', phantommask)
plt.imshow(phantommask)

xdim = 0.0459451172 # pixel size in object space 0.02297255 mm
ydim = 0.200 #pencilbeam step size in y direction mm



#attenuation map
PBattenuation_water = np.zeros((1024), dtype=np.float64)
PBattenuation_air = np.zeros((1024), dtype=np.float64)
PBenergy = 17.4 # 17.4 keV
mac_water = 1.258568 # cm^2/g for 17.4 keV x-ray photons
rho_water = 1.000E+00 #g/cm^3 density of water(phantom)
mac_air = 1.212672 # cm^2/g for 17.4 keV x-ray photons
rho_air = 1.205E-03 #g/cm^3 density of air




for i in xrange(1024):
    PBattenuation_water[i] = np.exp(-1*mac_water*rho_water*i*xdim*0.1)
    PBattenuation_air[i] = np.exp(-1*mac_air*rho_air*i*xdim*0.1)

#calculate [I_0 A t F_pix]
I0 = 1.31E+11 #Incident photon flux density Numbers/s/cm^2
A = 1.91E-04 #Beam area cm^2
t = 10 #dwell time 600 s

XrayFlux = np.zeros((250,1024), dtype=np.float64)
for i in xrange(250):
    if i< 85:
        XrayFlux[i,:] = I0 * A * t * PBattenuation_air[:]
    else:
        XrayFlux[i,:] = I0 * A * t * PBattenuation_water[:]

plt.imsave('XrayFlux.jpg', XrayFlux)
plt.imshow(XrayFlux)

#X-ray Fluorescence: Calculate [tau/rho_Y * p_k * rho_Y * deltax * omega_k * nu_alpha]
tau_Y = 96.41 # cm^2/g photoelectric cross-section
p_k = 0.8 # fraction of photoelectric interaction with K shell
rho_Y = 1 # g/cm^3 (estimate)
deltax = xdim*0.1 #voxel size 0.00459451172 cm
omega_k = 0.71 #fluorescence yield of k shell
nu_alpha = 0.854


density_y=np.zeros((250,2352), dtype=np.float64)
density_xy=np.zeros((250,1024), dtype=np.float64)

xrayfluorescence = np.zeros((250,1024), dtype=np.float64)
for i in xrange(250):
    for j in xrange(2352):
        
        density_y[i,j] = sum(tubemask[0+i*20:i*20+20, j])/20

      
for i in xrange(250):
    density_xy[i,0]=(density_y[i,0]+density_y[i,1]+density_y[i,2]*0.296875)/2.296875
    density_xy[i,1023]=(density_y[i,1023]+density_y[i,1022]+density_y[i,1021]*0.703125)/2.296875
    for j in xrange(1,1023):
        if j*2.296875-np.floor(j*2.296875) < 0.703125:
            density_xy[i,j] = (density_y[i,int(np.ceil(j*2.296875))]*(np.ceil(j*2.296875)-(j*2.296875)) + density_y[i,int(np.ceil(j*2.296875)+1)] + density_y[i,int(np.floor(2.296875*j+2.296875))]*(((j*2.296875+2.296875)-np.floor(2.296875*j+2.296875))))/2.296875
        else:    
            density_xy[i,j] = (density_y[i,int(np.ceil(j*2.296875))]*(np.ceil(j*2.296875)-(j*2.296875)) + density_y[i,int(np.ceil(j*2.296875)+1)] + density_y[i,int(np.ceil(j*2.296875)+2)] + density_y[i,int(np.floor(2.296875*j+2.296875))]*(((j*2.296875+2.296875)-np.floor(2.296875*j+2.296875))))/2.296875



for i in xrange(250):
    for j in xrange(1024):
        xrayfluorescence[i,j] = XrayFlux[i,j] * tau_Y * p_k * rho_Y * density_xy[i,j] * deltax * omega_k * nu_alpha


        
plt.imsave('xrayfluorescence.jpg', xrayfluorescence)
plt.imshow(xrayfluorescence)

print('Starting Detection of Xray Fluorescence')

#Detection of Xray Fluorescence
macfluoro_water = 1.74612 #cm^2/g
epsilon = 0.05 #detector quantum efficiency
geometric_efficiency = 0.00133106/1024/(4*np.pi)



path_x = np.zeros((1024,250), dtype=np.float64)
height_y = np.zeros((1024,250), dtype=np.float64)
for i in xrange(250):    
    if i > 84:
        d = (i-84)*ydim #mm
        for j in xrange(1024):
            if j < 512:
                h= (511-j)*xdim+xdim/2
                z=d*h/35.65
                path_x[j,i] = np.sqrt(z**2 + d**2) #mm
                height_y[j,i] = ((511-j)*0.027 + 0.0135)*d/56.60 #mm
            else:
                h= (j-512)*xdim+xdim/2
                z=d*h/35.65
                path_x[j,i] = np.sqrt(z**2 + d**2) #mm
                height_y[j,i] = ((j-512)*0.027 + 0.0135)*d/56.60 #mm

print('path_x and height_y calculated')                

path=np.zeros((1024,1024,165), dtype=np.float64)

for i in xrange(165):
    for j in xrange(1024):
        for k in xrange(1024):
            path[j,k,i]=0.1*np.sqrt(height_y[j,i+85]**2+path_x[k,i+85]**2) #cm
    print('finshed ' + str(i) + 'th itteration')
print('path calculated')


G_fluoro = np.exp(-1*macfluoro_water*rho_water*path)

Detector_XrayFluorescence =  np.zeros((1024,1024,165), dtype=np.float64)


print('calculateing Detector Fluorescence')

for i in xrange(165):
    for j in xrange(1024):
        for z in xrange(1024):
            Detector_XrayFluorescence[z,j,i]=xrayfluorescence[i+85,j] * G_fluoro[z,j,i] * epsilon * geometric_efficiency
    print('finshed ' + str(i) + 'th itteration')

plt.imsave('Dector_XrayFluorescence.jpg', Detector_XrayFluorescence[:,:,125])
plt.imshow(Detector_XrayFluorescence[:,:,125])

Detector_XrayFluorescence_samp1=Detector_XrayFluorescence[:,:,39]
Detector_XrayFluorescence_samp2=Detector_XrayFluorescence[:,:,40]
Detector_XrayFluorescence_samp3=Detector_XrayFluorescence[:,:,41]

print('calculating Random Signal')

RandomDetectorSignal=np.zeros((1024,1024,60), dtype=np.float64)
#Random number distribution
for z in xrange(60):
    for i in xrange(1024):
        for j in xrange(1024):
            RandomDetectorSignal[i,j,z]=np.random.poisson(Detector_XrayFluorescence_samp1[i,j]) + np.random.normal(0, 0.00063758)
    print('finshed ' + str(z) + 'th itteration of randomized detector noise')        
            

TotalRandDetectorSum=np.sum(RandomDetectorSignal, axis=2)
plt.imshow(TotalRandDetectorSum)
plt.imsave('TotalRandDetectorSum.jpg', TotalRandDetectorSum)

TotalRandDetectorMean=np.sum(RandomDetectorSignal, axis=2)
plt.imshow(TotalRandDetectorMean)
plt.imsave('TotalRandDetectorMean.jpg', TotalRandDetectorMean)

plt.imshow(RandomDetectorSignal[:,:,1])
plt.imsave('RandomDetectorSignal.jpg', RandomDetectorSignal[:,:,1])
