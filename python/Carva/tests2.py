#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Carva
"""

import numpy as np
import matplotlib.pyplot as plt

import sys, os
sys.path.insert(0, os.path.abspath('..'))
from tools.propagationTools import propagators
from tools.vortexTools import vortexProfiler

#%%----------
#System Parameters
#-----------------
# fundamentals
spatialSampling = 60.1e-3 # Divisor entero de pixSize (mm)
apertureRadius = 12.5 # Telescope - Lyot plane (mm)
wavelenght = 635e-6 # (mm)

# first lens (f_lens = 200 mm, r_lens = apertureRadius)
f_l1 = 200 # focal lenght (mm)
r_l1 = apertureRadius # lens radius (mm)

# SPP
shift = 0.5 # SPP axis position shift (mm)
pixSize = 60.1e-3 # SLM Pixel Pitch (mm)
Lvor = 2 # Topological Charge
NG = 2 # Number of gray levels

# second lens
f_l2 = f_l1 # lens 2 focal lenght (mm)
r_l2 = r_l1 # lens 2 radius (mm)
#%%-----------
#Tools
#-------------

VTools = vortexProfiler(dx=spatialSampling,radius=apertureRadius)
PTools = propagators(M = VTools.spaceSamples, dx =VTools.spatialStep, wl= wavelenght)
#%%-------------
#Entrance field
#---------------

U1 = VTools.placeAperture() # Uniform entrance field

#plt.imshow(np.abs(U1)) # -OK
#%%--------
#First lens
#----------

L1F = PTools.lens(f_l1) # Lens 1 function
U2 = U1*L1F # Post-lens_1 field

plt.imshow(np.angle(L1F)
#%%---------
#Prop to SPP
#-----------

A2 = PTools.FFT2D(U2)

plt.imshow(np.abs(A2))