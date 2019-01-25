#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Made by Grupo de Optica Aplicada - EAFIT
Juan José Cadavid Muñoz - 24-01-2019

Based on AFIC's work, this is an underdevelopment program which is refactorized
version into Object-Oriented Programming structure.
"""

import numpy as np
#import pyfftw
import matplotlib.pyplot as plt

#%% Class definition

class vortexProfiler:
    
    def __init__(self, dx = 26.e-3, p=12, radius=8.0):
        self.X = []
        self.Y = []
        self.x = []
        self.apertureRadius = radius #Apperture and Lyot stop radius
        self.spaceSamples = []
        self.halfSamples = []
        self.spaceSize = []
        self.halfSize = []
        self.factor = []
        self.spatialStep = dx #Modulator's pixel size in mm
        self.samplingPower = p
        
        self.computeSpace()
    
    def computeSpace(self):
        self.spaceSamples = 2**self.samplingPower #2**12=4096 #number of samples
        self.spaceSize = self.spaceSamples*self.spatialStep #100. #plane size in mm
        
        self.halfSamples = int(np.round(self.spaceSamples/2))
        self.halfSize = (self.spaceSize/2)
        
        self.x = np.arange(-self.halfSize,self.halfSize,self.spatialStep)
        self.X, self.Y = np.meshgrid(self.x,self.x) # Squared space is asumed
        self.apertureCenter = (apertureRadius/self.spatialStep)
        
        return self
            
    def createCircMask(self,a,b,step,w):
        x = np.arange(a,b-step,step); # Substract dx in order to obtain L instead of L+1
        X,Y = np.meshgrid(x,x)
        rho = np.sqrt(X**2 + Y**2)
        circMask = rho/abs(w) <= 0.5 # Bool circular shape descrption
        
        return circMask
    
    def placeAperture(self, centerPoint, windowCenter):
        
        window = np.zeros((self.spaceSamples, self.spaceSamples ))
        limM = int(centerPoint - windowCenter + 10) #Number of pixels needed to reach R2 plus a 
        limN = int(centerPoint + windowCenter + 10) # little extra-range, centered in M/2 (This process is made to avoid big calculations with the whole MxM matrix)
        window[limM:limN,limM:limN] = self.createCircMask(self.x[limM],self.x[limN],self.spatialStep,2*apertureRadius)
        
        return window
    
    def analyzeSpectrum(self,field):
        Ef = np.fft.fftshift(field)
        Ef = np.fft.fft2(Ef)
        Ef = np.fft.ifftshift(Ef)
        Ef = Ef*self.spatialStep**2
        
        return Ef
        
    def synthetizeSpectrum(self, field):
        Ei = np.fft.fftshift(field)
        Ei = np.fft.fft2(Ei)
        Ei = np.fft.ifftshift(Ei)
        Ei = Ei*1.0/self.spaceSize**2
        
        return Ei
    
    def cart2pol(self,x, y):
        rho = np.sqrt(x**2 + y**2)
        phi = np.arctan2(y, x)
        
        return(rho, phi)

    def SPP(self,Lvor,NG):
        
        x = np.arange(-self.halfSize, self.halfSize, self.spatialStep)
        xx, yy = np.meshgrid(x, x)
        phi = self.cart2pol(xx, yy)[1]
    
        phi1 = phi + np.pi #Matrix with entries within [0:2pi-step]
        phaseVor = np.mod(Lvor * phi1, 2*np.pi) #256-levels discretization
        
        phi2 = np.floor(phaseVor/(2*np.pi/NG)) #Matrix with whole-numbers between 0 and NG-1 
        phi3 = phi2/NG #phi3 is phi2 but normalized
        
        return np.exp(1j*(2*np.pi*phi3 - np.pi))

#%%-------------------------
#PROGRAM SETTINGS
#-------------------------
plotsEnabled = False

#%%---------------
#System Parameters
#-----------------
        
Lvor = 2 # Topologic Charge
NG = 10
spatialSampling = 26.e-3 # SLM Pixel Pitch (mm)
apertureRadius = 8.0 # Telescope - Lyot plane (mm)
#%%-------------------
#Vortex Analyzer Tools
#---------------------

vortexTools = vortexProfiler(dx=spatialSampling,radius=apertureRadius)
#%%----------------------
#Compute Field Properties
#------------------------

aperture = vortexTools.placeAperture(vortexTools.halfSamples,vortexTools.apertureCenter)
SLMfilterMask = vortexTools.SPP(Lvor,NG)
#%%--------------
#Propagate Fields
#----------------

SLMInput = vortexTools.analyzeSpectrum(aperture)
SLMPlane = SLMInput*SLMfilterMask
lyotPlane = vortexTools.synthetizeSpectrum(SLMPlane)
lyotAperturePlane = lyotPlane*aperture
outputField = vortexTools.analyzeSpectrum(lyotAperturePlane)
#%%---
#Plots
#-----

if plotsEnabled:

    viewRangeNA = vortexTools.halfSamples - 800
    viewRangeMA = vortexTools.halfSamples + 800
    
    viewRangeNB = vortexTools.halfSamples - 148
    viewRangeMB = vortexTools.halfSamples + 148
    
    viewRangeNC = vortexTools.halfSamples - 120
    viewRangeMC = vortexTools.halfSamples + 120
    
    # Truncated Input Plane at Stop
    plt.figure(2)
    ax=plt.axes()
    ax.set_title("$Telescope$ $Apperture,$ $z=0$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$")
    plt.imshow((abs(aperture)**2)[viewRangeNA:viewRangeMA,viewRangeNA:viewRangeMA],cmap='gray')
    #plt.show()
    
    # First Fourier Plane
    plt.figure(3)
    ax=plt.axes()
    ax.set_title("$First$ $Fourier$ $Plane,$ $z=f_{2}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$")
    plt.imshow((abs(SLMInput)**2)[viewRangeNB:viewRangeMB,viewRangeNB:viewRangeMB],cmap='gray')
    cb=plt.colorbar()
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Helicoidal Filter Mask
    plt.figure(4)
    ax=plt.axes()
    ax.set_title("$SPP.$ $m=%d,$ $N=%d,$ $z=f_{1}$"%(Lvor,NG),fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$")
    plt.imshow(np.angle(SLMfilterMask),cmap='gray')
    cb=plt.colorbar()
    cb.set_label(r"$Phase$ $value$",fontsize=12)
    #plt.show()
    
    # Signal with Filter at Modulating Plane
    fig=plt.figure(figsize=(14,8))
    
    ax0=plt.subplot(1,2,1)
    ax0.set_title("$Field$ $phase$ $after$ $SPP,$ $z=f_{1}$",fontsize=14,position=(0.5,1.0))
    ax0.set_xlabel("$N_{x}$",labelpad=8)
    ax0.set_ylabel("$N_{y}$") 
    map0=ax0.imshow(np.angle(SLMPlane)[viewRangeNC:viewRangeMC,viewRangeNC:viewRangeMC],cmap='gray')
    cb=plt.colorbar(map0,orientation='horizontal')
    cb.set_label(r"$Phase$ $value$",fontsize=12)
    
    ax1=plt.subplot(1,2,2)
    ax1.set_title("$Intensity$ $field$ $after$ $SPP,$ $z=f_{1}$",fontsize=14,position=(0.5,1.0))
    ax1.set_xlabel("$N_{x}$",labelpad=8)
    ax1.set_ylabel("$N_{y}$") 
    map1=ax1.imshow((abs(SLMPlane)**2)[viewRangeNC:viewRangeMC,viewRangeNC:viewRangeMC])
    cb=plt.colorbar(map1,orientation='horizontal')
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Lyot Plane - Filtered Signal - Inverse Fourier Field
    plt.figure(6)
    ax=plt.axes()
    ax.set_title("$Field$ $at$ $z=f_{1}+2*f_{2}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$") 
    plt.imshow((abs(lyotPlane)**2)[viewRangeNA:viewRangeMA,viewRangeNA:viewRangeMA],cmap='gray')
    cb=plt.colorbar()
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Lyod's Aperture - Signal truncation
    plt.figure(7)
    ax = plt.axes()
    ax.set_title("$Field$ $*$ $App.Stop,$ $z=f_{1}+2*f_{2}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$") 
    plt.imshow((abs(lyotAperturePlane)**2)[viewRangeNA:viewRangeMA,viewRangeNA:viewRangeMA],cmap='gray')
    cb=plt.colorbar(orientation='vertical')
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Focalized plane region - PSF
    
    plt.figure(8)
    ax = plt.axes()
    ax.set_title("$Intensity$ $field$ $at$ $camera,$ $z=f_{1}+2*f_{2}+2f_{3}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$") 
    plt.imshow((abs(outputField)**2)[viewRangeNC:viewRangeMC,viewRangeNC:viewRangeMC])#,cmap='gray')
    cb=plt.colorbar(orientation='vertical')
    cb.set_label(r"$Intensity$ $units$",fontsize=12)

    plt.show()