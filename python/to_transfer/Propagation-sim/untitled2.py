#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Made by Grupo de Optica Aplicada - EAFIT
Juan José Cadavid Muñoz - 24-01-2019

Based on AFIC's work, this is an underdevelopment program which is refactorized
version into Object-Oriented Programming structure.
"""

import numpy as np
import pyfftw
import matplotlib.pyplot as plt
from multiprocessing import cpu_count

#%% Class definition

class vortexProfiler:
    
    def __init__(self, dx = 26.0e-3, p=10, radius=8.0):
        self.X = []
        self.Y = []
        self.x = []
        self.rho = []
        self.phi = []
        self.apertureRadius = radius #Apperture and Lyot stop radius
        self.spaceSamples = []
        self.halfSamples = []
        self.spaceSize = []
        self.halfSize = []
        self.spatialStep = dx #Modulator's pixel size in mm
        self.samplingPower = p
        
        self.computeSpace() # Initialize space properties
    
    def computeSpace(self):
        self.spaceSamples = 2**self.samplingPower #2**12=4096 #number of samples
        self.spaceSize = self.spaceSamples*self.spatialStep #100. #plane size in mm
        self.halfSamples = int(np.fix(self.spaceSamples/2))
        self.halfSize = (self.spaceSize/2)
        
        self.x = np.arange(-self.halfSize,self.halfSize,self.spatialStep)
        self.X, self.Y = np.meshgrid(self.x,self.x) # Squared space is asumed
        self.rho = np.sqrt(self.X**2 + self.Y**2)
        self.phi = np.arctan2(self.Y, self.X)
        
        return self
            
    def createCircMask(self, w, relSize):
        circMask = self.rho/np.abs(w) <= relSize # Bool circular shape descrption
        return circMask
    
    def placeAperture(self):
        window = self.createCircMask(2*self.apertureRadius,0.5)
        return window
    
    def analyzeSpectrum(self,field):
        Ef = np.fft.fft2(field)
        Ef = Ef*self.spatialStep**2
        return Ef
        
    def synthetizeSpectrum(self, field):
        Ei = np.fft.ifft2(field)
        Ei = Ei*1.0/self.spaceSize**2
        return Ei

    def SPP(self,Lvor,NG):
        phi = Lvor*(self.phi + np.pi) #Matrix with entries within [0:2pi-step]
        phaseVor = np.mod(phi, 2*np.pi) #256-levels discretization
        
        phi = np.floor(phaseVor/(2*np.pi/NG)) #Matrix with whole-numbers between 0 and NG-1 
        phi = phi/NG #phi3 is phi2 but normalized
        return np.fft.fftshift(np.exp(1j*(2*np.pi*phi - np.pi)))
    
    def prepareFFTW(self,volumeSize):
        pyfftw.interfaces.cache.disable()

        fastLenghtT = pyfftw.next_fast_len(self.spaceSamples)
        fastLenghtV = pyfftw.next_fast_len(volumeSize)
        
        pyfftw.config.NUM_THREADS = cpu_count()
        pyfftw.config.PLANNER_EFFORT = 'FFTW_ESTIMATE'
        
        dataSize = (fastLenghtT,fastLenghtT,fastLenghtV)
        outputMatrixA = pyfftw.empty_aligned(dataSize, dtype='complex128', n=16)
        outputMatrixB = pyfftw.empty_aligned(dataSize, dtype='complex128', n=16)
        
        return outputMatrixA,outputMatrixB
    
    def plotData(self, viewRange, colHeader, rowHeader, dataSet, plotType, pad = 5):
        
        f, axes = plt.subplots(GLSize, TCSize, sharex=True, sharey=True)
        
        for ax, col in zip(axes[0], colHeader):
            ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                        xycoords='axes fraction', textcoords='offset points',
                        size='large', ha='center', va='baseline')
        
        for ax, row in zip(axes[:,0], rowHeader):
            ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                        xycoords=ax.yaxis.label, textcoords='offset points',
                        size='large', ha='right', va='center')
        
        if plotType in 'angle':
            for ax, index in zip(axes.flat, np.arange(0,volumeSize)):
                ax.imshow((np.angle(np.fft.fftshift(dataSet[:,:,index])[viewRange[0]:viewRange[1],viewRange[0]:viewRange[1]])),cmap='gray')
                
        elif plotType in 'intensity':
            for ax, index in zip(axes.flat, np.arange(0,volumeSize)):
                ax.imshow((np.abs(np.fft.fftshift(dataSet[:,:,index])[viewRange[0]:viewRange[1],viewRange[0]:viewRange[1]])**2),cmap='gray')
        
        f.tight_layout()
        f.subplots_adjust(left=0.15, top=0.95)
        plt.show()
        
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
NGStep = 2

spatialSampling = 60.1e-3 # SLM Pixel Pitch (mm)
apertureRadius = 4.0 # Telescope - Lyot plane (mm)
#%%-------------------
#Vortex Analyzer Tools
#---------------------

vortexTools = vortexProfiler(dx=spatialSampling,radius=apertureRadius)
#%%---------------
#Evaluation Ranges
#-----------------

TCRanges = np.arange(1,Lvor+1,TCStep)
GLRanges = np.arange(2,NG,2)
TCSize = len(TCRanges)
GLSize = len(GLRanges)
volumeSize = (TCSize)*(GLSize)

SLMPlanes = np.zeros((vortexTools.spaceSamples,vortexTools.spaceSamples,volumeSize),dtype = 'complex128')
allocatedMatrixSLM,allocatedMatrixOPT = vortexTools.prepareFFTW(volumeSize)
#%%----------------------
#Compute Field Properties
#------------------------

aperture = np.fft.fftshift(vortexTools.placeAperture())
SLMInput = vortexTools.analyzeSpectrum(aperture)

for grayIndex in range(0,GLSize):
    for TCIndex in range(0,TCSize):
        SLMPlanes[:,:,TCIndex + grayIndex*(TCSize)] = SLMInput*vortexTools.SPP(TCRanges[TCIndex],GLRanges[grayIndex])
#%%--------------
#Propagate Fields
#----------------
        
allocatedMatrixSLM[:] = SLMPlanes
allocatedMatrixOPT[:] = pyfftw.interfaces.numpy_fft.ifftn(allocatedMatrixSLM,axes=(0,1))*vortexTools.spatialStep**2
outputField = pyfftw.interfaces.numpy_fft.fftn(allocatedMatrixOPT,axes=(0,1))*1.0/vortexTools.spaceSize**2
#%%---
#Plots
#-----

if plotsEnabled:
    scaleRange = 0.25
    pixelShift = 1-(vortexTools.spaceSamples % 2) # Center graph if even matrix size is used
    viewRangeN = int((1 - scaleRange)*vortexTools.halfSamples) + pixelShift
    viewRangeM = int((1 + scaleRange)*vortexTools.halfSamples) + pixelShift
    
    cols = ['TC: {}'.format(col) for col in TCRanges]
    rows = ['NG: {}'.format(row) for row in GLRanges]
    
    vortexTools.plotData([viewRangeN,viewRangeM], cols, rows, allocatedMatrixOPT, 'intensity')
