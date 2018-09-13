"""
Binary mask for XF data
"""

import numpy as np

mask = np.load('XF_reconstruction/XF_recon.npy')
mask[mask < 0.9] = 0
mask[mask >= 0.9] = 1

np.save('XF_reconstruction/binarymaskXF.npy', mask)
