#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Carva
Last update in 26-03-2019
Coronografía UdeA-EAFIT

Propagation functions in object-oriented form
"""

import numpy as np

"""
Comments:

 - I think both classes could be merged into one, because they share atributes.
 
 - dx could be different from Modulator's pixel size, probably like half or quarter this size, to get better sampling in some cases for different functions.
"""
        
class propagators:

    """
    Planes are taken to be squares
    """
    
    def __init__(self, M, dx, wl):
        """
         - M: Number of samples in one side (Space samples) (M = L/dx, where L is side lenght).
                M, for example, could be len(u1), where u1 is the entrance field.
        
         - dx: source and observation plane side lenght delta.
         - wl: wavelenght.
        """
        
        # COMMON
        self.spaceSamples = M # Space samples
        self.halfSamples = M/2 # Hald space samples
        self.spaceSize = dx*M # Source and observation plane side lenght.
        self.halfSize = self.spaceSize/2 # Source and observation plane half side lenght.
        self.spatialStep = dx # Modulator's pixel size in mm ???    
        self.waveLenght = wl
        self.waveNumber = 2*np.pi/wl
        
        # lens
        self.R = []
        
        # TF prop and AS prop
        self.fx = []
        self.fX = []
        self.fY = []
        
        ## AS prop additional
        #self.fR = []
        
        # IR prop
        self.x = []
        self.X = []
        self.Y = []
        
        # FF prop
        self.L2 = []
        self.dx2 = []
        self.x2 = []
        self.X2 = []
        self.Y2 = []
        
    """
    The attributes out of "COMMON", are calculated only if necessary to optimize memory usage.
    Some of them are big matrixes so their effect wouldn't be neglectable.
    """
    
    def FFT2D(self, u): 
        """
        Does the 2D fft of a field, u, and returns it not-shifted
        """
        return np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(u)))
    
    def iFFT2D(self, A):
        """
        Does the inverse 2D fft of a frequency "field", A, and returns it not-shifted
        """
        return np.fft.ifftshift(np.fft.ifft2(np.fft.fftshift(A)))
    
    def CC(self):
        """
        Clear-caché function.
        
        Basically, cleans memory from not-indispensable variables.
        """
        self.R = []
        self.fx = []
        self.fX = []
        self.fY = []
        self.x = []
        self.X = []
        self.Y = []
        self.L2 = []
        self.dx2 = []
        self.x2 = []
        self.X2 = []
        self.Y2 = []
        
    def lens(self,f_lens):
        """
        REVIEW
        """
        if self.x == []:
            self.x =  np.arange(-self.spaceSize/2,self.spaceSize/2,self.spatialStep)
        if (self.X == []) or (self.Y == []):
            self.X, self.Y = np.meshgrid(self.x,self.x)
        if self.R == []:
            self.R = np.sqrt(self.X**2 + self.Y**2)
        
        Lop = np.exp(-(1j*np.pi)/(self.waveLenght*f_lens)*self.R**2)
        
        return Lop
    
    def TFOperator(self, z, cc = False):
        """
        Calculates TF propagation field-transformation function
         - z: propagation distance.
        """
        if self.fx == []:
            self.fx = np.arange(-1/(2*self.spatialStep),1/(2*self.spatialStep),1/self.spaceSize)
        if (self.fX == []) or (self.fY == []):
            self.fX,self.fY = np.meshgrid(self.fx,self.fx)
        
        TFop = np.exp(-1j*np.pi*self.waveLenght*z*(self.fX**2+self.fY**2)) # not-shifted operator
        
        if cc == True:
            self.CC()
        
        return TFop
        
    def IROperator(self, z, cc = False):
        """
        Calculates IR propagation field-transformation function
         - z: propagation distance.
        """
        if self.x == []:
            self.x =  np.arange(-self.spaceSize/2,self.spaceSize/2,self.spatialStep)
        if (self.X == []) or (self.Y == []):
            self.X, self.Y = np.meshgrid(self.x,self.x)
        
        IRop = np.fft.ifftshift(np.fft.fft2(np.fft.fftshift( 1/(1j*self.waveLenght*z)*np.exp(1j*self.waveNumber/(2*z)*(self.X**2+self.Y**2)) ))*self.spatialStep**2) # not-shifted operator
        
        if cc == True:
            self.CC()
        
        return IRop
        
    def FFProp(self, u, z, cc = False):
        """
        Calculates propagation to Fourier plane of u, complex optical field.
        
        returns u2, field at Fourier plane
        """
        
        if self.L2 == []:
            self.L2 = self.waveLenght*z/self.spatialStep
            self.dx2 = self.waveLenght*z/self.spaceSize
            self.x2 = np.arange(-self.L2/2,self.L2/2,self.dx2)
            
            self.X2,self.Y2 = np.meshgrid(self.x2,self.x2)
        
        c = 1/(1j*self.waveLenght*z)*np.exp(1j*self.waveNumber/(2*z)*(self.X2**2+self.Y2**2))
        u2 = self.FFT2D(u) * (c*self.spatialStep**2)
        
        if cc == True:
            self.CC()
            
        return u2

    def ASOperator(self,z, cc= False):
        """
        Calculates AS (angular espectrum) propagation field-transformation function
         - z: propagation distance
        """
        if self.fx == []:
            self.fx = np.arange(-1/(2*self.spatialStep),1/(2*self.spatialStep),1/self.spaceSize)
        if (self.fX == []) or (self.fY == []):
            self.fX, self.fY = np.meshgrid(self.fx,self.fx)
        
        ASop = np.exp(1j*self.waveNumber*z*np.sqrt(1-(self.waveLenght*self.fX)**2-(self.waveLenght*self.fY)**2))
        
        if cc == True:
            self.CC()
        
        return ASop
