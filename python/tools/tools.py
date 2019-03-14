#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Carva
Last update in 13-03-2019
Coronografía UdeA-EAFIT

Propagation functions in object-orientes form
"""

import numpy as np

"""
Comments:

 - I think both classes could be merged into one, because they share atributes. Thats why I put them in the same file (copying Juan José's work.
 
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
        self.spaceSize = dx/M # Source and observation plane side lenght.
        self.halfSize = self.spaceSize/2 # Source and observation plane half side lenght.
        self.spatialStep = dx # Modulator's pixel size in mm ???    
        self.waveLenght = wl
        self.waveNumber = 2*np.pi/wl
        
        # TF prop and AS prop
        self.fx = None
        self.fX = None
        self.fY = None
        
        ## AS prop additional
        #self.fR = None
        
        # IR prop
        self.x = None
        self.X = None
        self.Y = None
        
        # FF prop
        self.L2 = None
        self.dx2 = None
        self.x2 = None
        self.X2 = None
        self.Y2 = None
        
    """
    The attributes out of "COMMON", are calculated only if necessary to optimize memory usage. Some of them are big matrixes so their effect wouldn't be neglectable.
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
        self.fx = None
        self.fX = None
        self.fY = None
        self.x = None
        self.X = None
        self.Y = None
        self.L2 = None
        self.dx2 = None
        self.x2 = None
        self.X2 = None
        self.Y2 = None
        return None
    
    def TFOperator(self, z, cc = False):
        """
        Calculates TF propagation field-transformation function
         - z: propagation distance.
        """
        self.fx = np.arange(-1/(2*self.spatialStep),1/(2*self.spatialStep),1/self.spaceSize)
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
        self.x =  np.arange(-self.spaceSize/2,self.spaceSize/2,self.spatialStep)
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
        
        self.fx = np.arange(-1/(2*self.spatialStep),1/(2*self.spatialStep),1/self.spaceSize)
        
        self.fX, self.fY = np.meshgrid(self.fx,self.fx)
        
        ASop = np.exp(1j*self.waveNumber*z*np.sqrt(1-(self.waveLenght*self.fX)**2-(self.waveLenght*self.fY)**2))
        
        if cc == True:
            self.CC()
        
        return ASop