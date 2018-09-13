"""
Created on Thu Jan  4 14:01:35 2018

@author: bryanquigley
"""
from __future__ import division

import sys
sys.path.insert(0,"/home/bquigley/Research/XFXL Simulation 1.0.0")
import numpy as np
import matplotlib.image
import os

if not os.path.exists('XF_reconstruction'):
    os.makedirs('XF_reconstruction')

XF_recon = np.zeros((23,1024), dtype=np.float64)
XF_recon_80 = np.zeros((23,1024))
XF_recon_60 = np.zeros((23,1024))
XF_recon_40 = np.zeros((23,1024))
XF_recon_20 = np.zeros((23,1024))
XF_recon_10 = np.zeros((23,1024))
XF_recon_5  = np.zeros((23,1024))

for i in xrange(1,24):
    XFimage = np.load('XF_Simulated_images/Sim_XF_position_' + str(i) + '.npy')
    XF_recon[i-1,:]=np.sum(XFimage, axis = 0)
    XF_recon_80[i-1,:]=np.sum(XFimage[102:921,:], axis = 0)
    XF_recon_60[i-1,:]=np.sum(XFimage[204:819,:], axis = 0)
    XF_recon_40[i-1,:]=np.sum(XFimage[307:717,:], axis = 0)
    XF_recon_20[i-1,:]=np.sum(XFimage[410:614,:], axis = 0)
    XF_recon_10[i-1,:]=np.sum(XFimage[461:563,:], axis = 0)
    XF_recon_5[i-1,:] =np.sum(XFimage[486:538,:], axis = 0)


np.save('XF_reconstruction/XF_sim_recon_100.npy', XF_recon)
np.save('XF_reconstruction/XF_sim_recon_80.npy', XF_recon_80)
np.save('XF_reconstruction/XF_sim_recon_60.npy', XF_recon_60)
np.save('XF_reconstruction/XF_sim_recon_40.npy', XF_recon_40)
np.save('XF_reconstruction/XF_sim_recon_20.npy', XF_recon_20)
np.save('XF_reconstruction/XF_sim_recon_10.npy', XF_recon_10)
np.save('XF_reconstruction/XF_sim_recon_05.npy', XF_recon_5)

np.save('XF_reconstruction/XF_recon', XF_recon)
matplotlib.image.imsave('XF_reconstruction/XF_recon.tif', np.array(XF_recon), cmap = 'gray')

np.save('XF_reconstruction/XF_recon_80', XF_recon_80)
matplotlib.image.imsave('XF_reconstruction/XF_recon_80.tif', np.array(XF_recon_80), cmap = 'gray')

np.save('XF_reconstruction/XF_recon_60', XF_recon_60)
matplotlib.image.imsave('XF_reconstruction/XF_recon_60.tif', np.array(XF_recon_60), cmap = 'gray')

np.save('XF_reconstruction/XF_recon_40', XF_recon_40)
matplotlib.image.imsave('XF_reconstruction/XF_recon_40.tif', np.array(XF_recon_40), cmap = 'gray')

np.save('XF_reconstruction/XF_recon_20', XF_recon_20)
matplotlib.image.imsave('XF_reconstruction/XF_recon_20.tif', np.array(XF_recon_20), cmap = 'gray')

np.save('XF_reconstruction/XF_recon_10', XF_recon_10)
matplotlib.image.imsave('XF_reconstruction/XF_recon_10.tif', np.array(XF_recon_10), cmap = 'gray')
