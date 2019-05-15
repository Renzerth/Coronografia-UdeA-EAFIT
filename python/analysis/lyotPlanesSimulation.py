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

import numpy as np
from scipy.io import savemat
        
#%%--------------
#PROGRAM SETTINGS
#----------------
        
plotsEnabled = True
#%%---------------
#System Parameters
#-----------------
        
Lvor = 10 # Topologic Charge
TCStep = 1
NGmin = 2
NGmax = 10
NG = 4

spatialSampling = 20.1e-3 # SLM Pixel Pitch (mm)
apertureRadius = 2.0 # Telescope - Lyot plane (mm)
#%%-------------------
#Vortex Analyzer Tools
#---------------------

vortexTools = vortexProfiler(dx=spatialSampling,p=10,radius=apertureRadius)
#%%---------------
#Evaluation Ranges
#-----------------

TCRanges = np.arange(1,Lvor+1,TCStep)

GLRanges = np.fix(np.linspace(NGmin,NGmax,NG)).astype('int')
#GLRanges = np.array([12,16,24,32,64,128,256])

TCSize = len(TCRanges)
GLSize = len(GLRanges)
volumeSize = (TCSize)*(GLSize)

SLMPlanes = np.zeros((vortexTools.spaceSamples,vortexTools.spaceSamples,volumeSize),dtype = 'complex128')
SLMPMasks = np.zeros((vortexTools.spaceSamples,vortexTools.spaceSamples,volumeSize),dtype = 'complex128')
allocatedMatrixSLM,allocatedMatrixLyot = vortexTools.prepareFFTW(volumeSize)
#%%----------------------
#Compute Field Properties
#------------------------

aperture = np.fft.fftshift(vortexTools.placeAperture(0.5))
#aperture = np.fft.fftshift(vortexTools.placeGaussianAperture(1.0))

lyotAperture = np.fft.fftshift(vortexTools.placeAperture(0.5))
SLMInput = vortexTools.analyzeSpectrum(aperture)

for TCIndex in range(0,TCSize):
    for grayIndex in range(0,GLSize):

        singleIndex = TCIndex + grayIndex*(TCSize)
        SLMPMasks[:,:,singleIndex] = vortexTools.SPP(TCRanges[TCIndex],GLRanges[grayIndex])
        SLMPlanes[:,:,singleIndex] = SLMInput*SLMPMasks[:,:,singleIndex]
#%%----------------------
#Generate PSF Reference
#------------------------
PSFreference = vortexTools.analyzeSpectrum(vortexTools.synthetizeSpectrum(SLMInput)*lyotAperture);
#%%--------------
#Propagate Fields
#----------------
        
allocatedMatrixSLM[:] = SLMPlanes
allocatedMatrixLyot[:] = vortexTools.propagateField(allocatedMatrixSLM,'forward') # Lyot's Plane Response
allocatedMatrixLyotTrunc = allocatedMatrixLyot*lyotAperture[:,:,np.newaxis]
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
    vortexTools.plotData([viewRangeN,viewRangeM], cols, rows, allocatedMatrixLyotTrunc, 'intensity')
    vortexTools.plotData([viewRangeN,viewRangeM], cols, rows, PSFoutputFields, 'log')
#%%-------------
#Save Data Files
#---------------
matlaborpyton = 1 # 1: MATLAB; 2: Python
mdict={'PSFoutputFields':PSFoutputFields,'PSFreference':PSFreference,'TCRanges':TCRanges,'GLRanges':GLRanges}
saveDir = '../data/Simulations/PSFs'

if matlaborpyton == 1:
    savemat(saveDir + '.mat',mdict)
else:
    np.save(saveDir + '.npy',mdict)