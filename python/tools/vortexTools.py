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
        self.apertureRadius = radius # Apperture and Lyot stop radius
        self.spaceSamples = []
        self.halfSamples = []
        self.spaceSize = []
        self.halfSize = []
        self.spatialStep = dx # Modulator's pixel size in mm
        self.samplingPower = p
        
        self.computeSpace() # Initialize space properties
    
    def computeSpace(self):
        self.spaceSamples = 2**self.samplingPower #2**12=4096 #number of samples
        self.spaceSize = self.spaceSamples*self.spatialStep #100. # Plane size in mm
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
    
    def placeAperture(self,radius=0.5):
        window = self.createCircMask(2*self.apertureRadius,radius).astype('double')
        return window
    
    def placeGaussianAperture(self,sigma=0.5):
        window = np.exp(-(self.rho/sigma)**2)
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
        phi = Lvor*(self.phi + np.pi) # Matrix with entries within [0:2pi-step]
        phaseVor = np.mod(phi, 2*np.pi) # 256-levels discretization
        
        phi = np.floor(phaseVor/(2*np.pi/NG)) # Matrix with whole-numbers between 0 and NG-1 
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
    
    def propagateField(self, fields, direction):
        if direction in 'forward':
            fieldTransformation = pyfftw.interfaces.numpy_fft.fftn(fields,axes=(0,1))*self.spatialStep**2
        elif direction in 'backward':
            fieldTransformation= pyfftw.interfaces.numpy_fft.ifftn(fields,axes=(0,1))*1.0/self.spaceSize**2
            
        return fieldTransformation
    
    def plotData(self, viewRange, colHeader, rowHeader, dataSet, plotType, pad = 5):
        
        f, axes = plt.subplots(len(rowHeader), len(colHeader), sharex=True, sharey=True)
        
        for ax, col in zip(axes[0], colHeader):
            ax.annotate(col, xy=(0.5, 1), xytext=(0, pad),
                        xycoords='axes fraction', textcoords='offset points',
                        size=20, ha='center', va='baseline')
        
        for ax, row in zip(axes[:,0], rowHeader):
            ax.annotate(row, xy=(0, 0.5), xytext=(-ax.yaxis.labelpad - pad, 0),
                        xycoords=ax.yaxis.label, textcoords='offset points',
                        size=20, ha='right', va='center',rotation='vertical')
        
        if plotType in 'angle':
            for index, ax in enumerate(axes.flat,0):
                ax.imshow((np.angle(np.fft.fftshift(dataSet[:,:,index])[viewRange[0]:viewRange[1],viewRange[0]:viewRange[1]])),cmap='gray')
                
        elif plotType in 'intensity':
            for index, ax in enumerate(axes.flat,0):
                ax.imshow((np.abs(np.fft.fftshift(dataSet[:,:,index])[viewRange[0]:viewRange[1],viewRange[0]:viewRange[1]])**2),cmap='gray')
        elif plotType in 'log':
            for index, ax in enumerate(axes.flat,0):
                ax.imshow((np.log10(np.abs(np.fft.fftshift(dataSet[:,:,index])[viewRange[0]:viewRange[1],viewRange[0]:viewRange[1]])**2)))
        
        f.tight_layout()
        f.subplots_adjust(left=0.15, top=0.95)
        plt.show()