#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Made by Grupo de Optica Aplicada - EAFIT
Juan José Cadavid Muñoz - 24-01-2019

Based on AFIC's work, this is an underdevelopment program which is refactorized
version into Object-Oriented Programming structure.
"""

import sys, os; 
sys.path.insert(0, os.path.abspath('..'))
from tools.vortexTools import vortexProfiler
from scipy import tensordot

import numpy as np
        
#%%--------------
#PROGRAM SETTINGS
#----------------
        
plotsEnabled = True
#%%---------------
#System Parameters
#-----------------
        
Lvor = 6 # Topologic Charge
TCStep = 1
NG = 10
NGStep = 1

spatialSampling = 60.1e-3 # SLM Pixel Pitch (mm)
apertureRadius = 4.0 # Telescope - Lyot plane (mm)
#%%-------------------
#Vortex Analyzer Tools
#---------------------

vortexTools = vortexProfiler(dx=spatialSampling,p=9,radius=apertureRadius)
#%%---------------
#Evaluation Ranges
#-----------------

TCRanges = np.arange(1,Lvor+1,TCStep)
GLRanges = np.arange(2,NG,2)
TCSize = len(TCRanges)
GLSize = len(GLRanges)
volumeSize = (TCSize)*(GLSize)

SLMPlanes = np.zeros((vortexTools.spaceSamples,vortexTools.spaceSamples,volumeSize),dtype = 'complex128')
SLMPMasks = np.zeros((vortexTools.spaceSamples,vortexTools.spaceSamples,volumeSize),dtype = 'complex128')
allocatedMatrixSLM,allocatedMatrixLyot = vortexTools.prepareFFTW(volumeSize)
#%%----------------------
#Compute Field Properties
#------------------------

aperture = np.fft.fftshift(vortexTools.placeAperture())
SLMInput = vortexTools.analyzeSpectrum(aperture)
lyotAperture = aperture[:,:,np.newaxis]

for grayIndex in range(0,GLSize):
    for TCIndex in range(0,TCSize):
        SLMPMasks[:,:,TCIndex + grayIndex*(TCSize)] = vortexTools.SPP(TCRanges[TCIndex],GLRanges[grayIndex])
        SLMPlanes[:,:,TCIndex + grayIndex*(TCSize)] = SLMInput*SLMPMasks[:,:,TCIndex + grayIndex*(TCSize)]
#%%--------------
#Propagate Fields
#----------------
        
allocatedMatrixSLM[:] = SLMPlanes
allocatedMatrixLyot[:] = vortexTools.propagateField(allocatedMatrixSLM,'forward') # Lyot's Plane Response
allocatedMatrixLyotTrunc = allocatedMatrixLyot*aperture[:,:,np.newaxis]
PSFoutputFields = vortexTools.propagateField(allocatedMatrixLyotTrunc,'backward') # PSF response
#%%---
#Plots
#-----

if plotsEnabled:
    scaleRange = 0.35 # 1 for no zoom -> 0 for one pixel zoom
    pixelShift = 1-(vortexTools.spaceSamples % 2) # Center graph if even matrix size is used
    viewRangeN = int((1 - scaleRange)*vortexTools.halfSamples) + pixelShift
    viewRangeM = int((1 + scaleRange)*vortexTools.halfSamples) + pixelShift
    
    cols = ['TC: {}'.format(col) for col in TCRanges]
    rows = ['NG: {}'.format(row) for row in GLRanges]
    
    vortexTools.plotData([viewRangeN,viewRangeM], cols, rows, SLMPMasks, 'angle')
    vortexTools.plotData([viewRangeN,viewRangeM], cols, rows, allocatedMatrixLyot, 'intensity')
    vortexTools.plotData([viewRangeN,viewRangeM], cols, rows, PSFoutputFields, 'log')
