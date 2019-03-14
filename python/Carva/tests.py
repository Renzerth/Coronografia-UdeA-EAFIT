#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Carva
"""

import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.insert(0, os.path.abspath('..'))
from tools.tools import propagationTools

#%%------------
# Input field
#--------------

L = 4 # cm
M = 1024
dx = L/M
wl = 0.633e-6 # cm

x, y = np.arange(-L/2,L/2,dx), np.arange(-L/2,L/2,dx)
X,Y = np.meshgrid(x,y)


a = 15 # some const. for the field

U1 = np.exp(-a*(X**2+Y**2))

# Now if we want a circular stop or radiud r = 0.5, is as simple as:

LS = np.sqrt(X**2+Y**2) < 1.5 # boolean-like array for simplifying the calculations

# Entry stop
U1 = U1*LS # U1 is the field at entry stop

del L, x, y, X, Y, a

#%%-------------
# test on propagators
#---------------

propTools = propagationTools(M, dx, wl)

A1 = propTools.iFFT2D(U1)
PF = propTools.TFOperator(900000,cc=False)

A2 = A1*PF
U2 = propTools.iFFT2D(A2)

#print(propTools.spaceSize)
#plt.imshow(np.abs(U1))
plt.imshow(np.abs(U2))
