# -*- coding: utf-8 -*-
"""
Created on Wed Nov 10 23:33:33 2021

@author: arefk
"""


from __future__ import print_function  # print('me') instead of print 'me'
from __future__ import division  # 1/2 == 0.5, not 0


# - import common modules
import numpy as np  # the Python array package


def mutualInfo(Im1,Im2):

    t1_slice = Im1
    t2_slice = Im2

    hist_2d, x_edges, y_edges = np.histogram2d(t1_slice.ravel(),t2_slice.ravel(),bins=20)

    hist_2d_log = np.zeros(hist_2d.shape)
    non_zeros = hist_2d != 0
    hist_2d_log[non_zeros] = np.log(hist_2d[non_zeros])
    
    pxy = hist_2d / float(np.sum(hist_2d))
    px = np.sum(pxy, axis=1) # marginal for x over y
    py = np.sum(pxy, axis=0) # marginal for y over x
    px_py = px[:, None] * py[None, :] # Broadcast to multiply marginals
    # Now we can do the calculation using the pxy, px_py 2D arrays
    nzs = pxy > 0 # Only non-zero pxy values contribute to the sum
    MI = np.sum(pxy[nzs] * np.log(pxy[nzs] / px_py[nzs]))
    
    return MI
