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
<<<<<<< HEAD
wavelength = 635e-6 # (mm)
propc = 0

# first lens (f_lens = 200 mm, r_lens = apertureRadius)
f_l1 = 200 # focal length (mm)
=======
wavelenght = 635e-6 # (mm)
propc = 0

# first lens (f_lens = 200 mm, r_lens = apertureRadius)
f_l1 = 200 # focal lenght (mm)
>>>>>>> origin/master
r_l1 = apertureRadius # lens radius (mm)

# SPP
shift = 35.0 # SPP axis position shift (mm)
pixSize = 60.1e-3 # SLM Pixel Pitch (mm)
Lvor = 2 # Topological Charge
NG = 2 # Number of gray levels

# second lens
<<<<<<< HEAD
f_l2 = f_l1 # lens 2 focal length (mm)
=======
f_l2 = f_l1 # lens 2 focal lenght (mm)
>>>>>>> origin/master
r_l2 = r_l1 # lens 2 radius (mm)


# first prop
distanceOne = f_l1-shift

# second prop
distanceTwo = f_l1+shift+f_l2

# thrid prop
distanceThree = f_l2

#%%-----------
#Tools
#-------------

VTools = vortexProfiler(dx=spatialSampling,radius=apertureRadius, p = 12)
<<<<<<< HEAD
PTools = propagators(M = VTools.spaceSamples, dx =VTools.spatialStep, wl= wavelength)
=======
PTools = propagators(M = VTools.spaceSamples, dx =VTools.spatialStep, wl= wavelenght)
>>>>>>> origin/master
#%%-------------------------------------------------
#Propagation (Memory responsible condensed ver.)
#---------------------------------------------------

U1 = VTools.placeAperture() # Uniform entrance field

L1F = PTools.lens(f_l1) # Lens 1 function

U1 = U1*L1F # Post-lens_1 field

del L1F

A1 = PTools.FFT2D(U1) # U2 to frequency space

del U1

if propc == 0:
    POperator1 = PTools.TFOperator(distanceOne)
elif propc == 1:
    POperator1 = PTools.ASOperator(distanceOne)
elif propc == 2:
    POperator1 =PTools.IROperator(distanceOne)

A1 = A1*POperator1

del POperator1

U1 = PTools.iFFT2D(A1)

del A1

SPP = np.fft.ifftshift(VTools.SPP(Lvor, NG))

U1 = U1*SPP
A1 = PTools.FFT2D(U1)

del U1

if propc == 0:
    POperator2 = PTools.TFOperator(distanceTwo)
elif propc == 1:
    POperator2 = PTools.ASOperator(distanceTwo)
elif propc == 2:
    POperator2 =PTools.IROperator(distanceTwo)
    
A1 = A1*POperator2

del POperator2

U1 = PTools.iFFT2D(A1)

del A1

L2F = PTools.lens(f_l2)
U1 = U1*L2F

A1 = PTools.FFT2D(U1)

del U1

if propc == 0:
    POperator3 = PTools.TFOperator(distanceThree)
elif propc == 1:
    POperator3 = PTools.ASOperator(distanceThree)
elif propc == 2:
    POperator3 =PTools.IROperator(distanceThree)


A1 = A1*POperator3

del POperator3

U1 = PTools.iFFT2D(A1)

del A1

#%%

plt.imshow(np.angle(U1),'gray')
plt.colorbar()

#plt.savefig('methodIR.png')

#%%----------------------------------------------------------------------------
###############################################################################
############################## HISTORY ########################################
###############################################################################

"""
#%%-------------
#Entrance field
#---------------

U1 = VTools.placeAperture() # Uniform entrance field

plt.imshow(np.abs(U1)) # -OK

#%%
plt.imshow(np.angle(U1))
#%%
<<<<<<< HEAD
PTools.wavelength
=======
PTools.waveLenght
>>>>>>> origin/master
#%%--------
#First lens
#----------

L1F = PTools.lens(f_l1) # Lens 1 function
U2 = U1*L1F # Post-lens_1 field

plt.imshow(np.abs(U2))
plt.colorbar()
#%%
plt.imshow(np.angle(L1F))
plt.colorbar()
#%%---------
#Prop to SPP
#-----------

A2 = PTools.FFT2D(U2) # U2 to frequency space

# which propagator to use?
# in this version, im not considering optimization.

propc = 2
distanceOne = f_l1-shift

if propc == 0:
    POperator1 = PTools.TFOperator(distanceOne)
elif propc == 1:
    POperator1 = PTools.ASOperator(distanceOne)
elif propc == 2:
    POperator1 =PTools.IROperator(distanceOne)

A3 = A2*POperator1

U3 = PTools.iFFT2D(A3)

plt.imshow(np.abs(U3))
#%% phase plot
plt.imshow(np.angle(U3))
#%%-
#SPP
#---

SPP = np.fft.ifftshift(VTools.SPP(Lvor, NG))

U4 = U3*SPP
A4 = PTools.FFT2D(U4)

plt.imshow(np.abs(U4))
#%% SPP plot

plt.imshow(np.angle(SPP),'gray')
#%% phase plot
plt.imshow(np.angle(U4))
#%%----------------
#Second propagation
#------------------

distanceTwo = f_l1+shift+f_l2


if propc == 0:
    POperator2 = PTools.TFOperator(distanceTwo)
elif propc == 1:
    POperator2 = PTools.ASOperator(distanceTwo)
elif propc == 2:
    POperator2 =PTools.IROperator(distanceTwo)
    
A5 = A4*POperator2
U5 = PTools.iFFT2D(A5)

plt.imshow(np.abs(U5))
#%% phase plot
plt.imshow(np.angle(U5))
#%%---------
#Second lens
#-----------

L2F = PTools.lens(f_l2)
U6 = U5*L2F

plt.imshow(np.abs(U6))
#%% phase plot
plt.imshow(np.angle(U6))
#%%---------------
#Final propagation
#-----------------

distanceThree = f_l2

if propc == 0:
    POperator3 = PTools.TFOperator(distanceThree)
elif propc == 1:
    POperator3 = PTools.ASOperator(distanceThree)
elif propc == 2:
    POperator3 =PTools.IROperator(distanceThree)

A6 = PTools.FFT2D(U6)

A7 = A6*POperator3
U7 = PTools.iFFT2D(A7)
#%% final plots
plt.imshow(np.abs(U7))
#%% phase plot
plt.imshow(np.angle(U7))

#%%---------------
#Condensed program V1
#-----------------

U1 = VTools.placeAperture() # Uniform entrance field

L1F = PTools.lens(f_l1) # Lens 1 function

U2 = U1*L1F # Post-lens_1 field
A2 = PTools.FFT2D(U2) # U2 to frequency space

if propc == 0:
    POperator1 = PTools.TFOperator(distanceOne)
elif propc == 1:
    POperator1 = PTools.ASOperator(distanceOne)
elif propc == 2:
    POperator1 =PTools.IROperator(distanceOne)

A3 = A2*POperator1

U3 = PTools.iFFT2D(A3)

SPP = np.fft.ifftshift(VTools.SPP(Lvor, NG))

U4 = U3*SPP
A4 = PTools.FFT2D(U4)

if propc == 0:
    POperator2 = PTools.TFOperator(distanceTwo)
elif propc == 1:
    POperator2 = PTools.ASOperator(distanceTwo)
elif propc == 2:
    POperator2 =PTools.IROperator(distanceTwo)
    
A5 = A4*POperator2
U5 = PTools.iFFT2D(A5)

L2F = PTools.lens(f_l2)
U6 = U5*L2F

if propc == 0:
    POperator3 = PTools.TFOperator(distanceThree)
elif propc == 1:
    POperator3 = PTools.ASOperator(distanceThree)
elif propc == 2:
    POperator3 =PTools.IROperator(distanceThree)

A6 = PTools.FFT2D(U6)

A7 = A6*POperator3
U7 = PTools.iFFT2D(A7)
"""