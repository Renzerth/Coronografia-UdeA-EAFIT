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
    
    def __init__(self, spaceSize, spaceSamples):
        self.X = []
        self.Y = []
        self.x = []
        self.rho = []
        self.phi = []
        self.phiB = []
        self.spaceSamples = spaceSamples
        self.halfSamples = []
        self.spaceSize = spaceSize
        self.halfSize = []
        self.spatialStep = spaceSize/spaceSamples # Space's pixel pitch in mm
        
        self.computeSpace() # Initialize space properties
    
    def computeSpace(self):
        self.halfSamples = int(np.floor((self.spaceSamples+1)/2))
        self.halfSize = (self.spaceSize/2)
        
        self.x = np.arange(-self.halfSize,self.halfSize,self.spatialStep)
        self.X, self.Y = np.meshgrid(self.x,self.x) # Squared space is asumed
        self.rho = np.sqrt(self.X**2 + self.Y**2)
        self.phiB = np.arctan2(self.Y, self.X)
        self.generateCSpiral()
        
        return self
            
    def createCircMask(self, circWidth, relativeSize):
        circMask = self.rho/np.abs(circWidth) < relativeSize # Bool circular shape descrption
        return circMask
    
    def placeAperture(self, apertureDiameter, relativeSize=1):
        window = self.createCircMask(apertureDiameter,relativeSize).astype('double')
        return window
    
    def placeGaussianAperture(self,sigma=0.5):
        window = np.exp(-(self.rho/sigma)**2)
        return window
    
    def generateCSpiral(self):
        m = self.spaceSamples # Row Size
        n = self.spaceSamples # Column Size
        
        Row = np.arange(0,m,1) - m/2 # All Rows
        
        Column1 = np.arange(0,n/2+1,1) - n/2 # First left half
        Column2 = np.arange(n/2,n,1) - n/2 # Last right half
        
        x1,y1 = np.meshgrid(Column1,Row)
        th1 = np.arctan(y1/x1) + np.pi/2 # First Half Angular transition

        x2,y2 = np.meshgrid(Column2,Row)
        th2 = np.arctan(y2/x2) + np.pi + np.pi/2 # Last Half Angular transition
        
        phaseMask = np.exp(1j*np.concatenate((th1[:,0:(n//2)],th2[:,0:n-1]),axis=1))
        phaseMask[self.halfSamples,self.halfSamples] = -np.pi
        self.phi = np.angle(phaseMask)
        return self
    
    def analyzeSpectrum(self,field):
        Ef = np.fft.fft2(field)
        Ef = Ef*self.spatialStep**2
        return Ef
        
    def synthetizeSpectrum(self, field):
        Ei = np.fft.ifft2(field)
        Ei = Ei*1.0/self.spaceSize**2
        return Ei

    def discretizeSPP(self,TC,NG):
        phi = TC*(self.phi + np.pi) # Matrix with entries within [0:2pi-step]
        phaseVor = np.mod(phi, 2*np.pi) # 256-levels discretization
        
        phi = np.floor(phaseVor/(2*np.pi/NG)) # Matrix with whole-numbers between 0 and NG-1 
        phi = phi/NG #phi is normalized
        vortexMask = np.exp(1j*(2*np.pi*phi - np.pi))
        
        return np.fft.fftshift(vortexMask)
    
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
    
    def propagateField(self, fields, direction):
        if direction in 'forward':
            fieldTransformation = pyfftw.interfaces.numpy_fft.fftn(fields,axes=(0,1))*self.spatialStep**2
        elif direction in 'backward':
            fieldTransformation= pyfftw.interfaces.numpy_fft.ifftn(fields,axes=(0,1))*1.0/self.spaceSize**2
            
        return fieldTransformation
    
    def plotData(self, viewRange, colHeader, rowHeader, dataSet, plotType, pad = 5):
        
        dataExtend = [-self.halfSize,self.halfSize,-self.halfSize,self.halfSize];
        f, axes = plt.subplots(len(rowHeader), len(colHeader), sharex=True, sharey=True)
        
        for ax, col in zip(axes[0], colHeader):
            ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                        xycoords='axes fraction', textcoords='offset points',
                        size=20, ha='center', va='baseline')
        
        for ax, row in zip(axes[:,0], rowHeader):
            ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                        xycoords=ax.yaxis.label, textcoords='offset points',
                        size=20, ha='right', va='center', rotation='vertical')
        
        if plotType in 'angle':
            for index, ax in enumerate(axes.flat,0):
                ax.imshow((np.angle(np.fft.fftshift(dataSet[:,:,index])[viewRange[0]:viewRange[1],viewRange[0]:viewRange[1]])), cmap='gray', interpolation='none', extent=dataExtend)
                ax.set_aspect('equal', 'box')
                
        elif plotType in 'intensity':
            for index, ax in enumerate(axes.flat,0):
                ax.imshow((np.abs(np.fft.fftshift(dataSet[:,:,index])[viewRange[0]:viewRange[1],viewRange[0]:viewRange[1]])**2),cmap='gray', interpolation='none', extent=dataExtend)
                ax.set_aspect('equal', 'box')
                
        elif plotType in 'log':
            for index, ax in enumerate(axes.flat,0):
                ax.imshow((np.log10(np.abs(np.fft.fftshift(dataSet[:,:,index])[viewRange[0]:viewRange[1],viewRange[0]:viewRange[1]])**2)), interpolation='none', extent=dataExtend)
                ax.set_aspect('equal', 'box')
        
        f.tight_layout()
        f.subplots_adjust(left=0.15, top=0.95)
        plt.show()