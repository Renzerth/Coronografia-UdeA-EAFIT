#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Made by Grupo de Optica Aplicada - EAFIT
Juan José Cadavid Muñoz - 24-01-2019

Based on AFIC's work, this is an underdevelopment program which is refactorized
version into Object-Oriented Programming structure.
"""

import numpy as np
import matplotlib.pyplot as plt

import sys, os; 
sys.path.insert(0, os.path.abspath('..'))
from tools.vortexTools import vortexProfiler
np.seterr(divide='ignore', invalid='ignore') # avoid np.arctan(0/0) notification

#%%--------------
#PROGRAM SETTINGS
#----------------
        
plotsEnabled = True
#%%---------------
#System Parameters
#-----------------
        
Lvor = 1 # Topological Charge
NG = 256 # Number of gray levels [1,256]
spatialSampling = 120.1e-3 # Space sampling size (mm)
apertureRadius = 1.0 # Telescope - Lyot plane (mm)
LyotApertureRadius = 1.0
#%%-------------------
#Vortex Analyzer Tools
#---------------------

vortexTools = vortexProfiler(dx=spatialSampling,p=10,radius=apertureRadius)
#%%----------------------
#Compute Field Properties
#------------------------

telescopeAperture = vortexTools.placeAperture();
lyotAperture = vortexTools.placeAperture(LyotApertureRadius)
SLMfilterMask = vortexTools.discretizeSPP(Lvor,NG)
#SLMfilterMask = np.fft.fftshift(np.exp(1j*vortexTools.phiB*Lvor))
#%%--------------
#Propagate Fields
#----------------

SLMInput = vortexTools.analyzeSpectrum(telescopeAperture)
SLMPlane = SLMInput*SLMfilterMask
lyotPlane = vortexTools.synthetizeSpectrum(SLMPlane)
lyotAperturePlane = lyotPlane*lyotAperture
outputField = vortexTools.analyzeSpectrum(lyotAperturePlane)
#%%---
#Plots
#-----

if plotsEnabled:
    
    spatialCoords = vortexTools.x
    dx = spatialSampling
    spatialExtent = [spatialCoords[0]-dx, spatialCoords[-1]+dx, spatialCoords[0]-dx, spatialCoords[-1]+dx]
    
    scaleRange = 1 # 1 -> To full field, 0.1 -> To close up, values < 0.1 renders nothing
    pixelShift = 1-(vortexTools.spaceSamples % 2) # Center graph if even matrix size is used
    
    viewRangeN = int((1 - scaleRange)*vortexTools.halfSamples)+ pixelShift
    viewRangeM = int((1 + scaleRange)*vortexTools.halfSamples) + pixelShift
    
    # Truncated Input Plane at Stop
    plt.figure(2)
    ax=plt.axes()
    ax.set_title("$Telescope$ $Apperture,$ $z=0$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$ [mm]",labelpad=8)
    ax.set_ylabel("$N_{y}$ [mm]")
    plt.imshow(((telescopeAperture[viewRangeN:viewRangeM,viewRangeN:viewRangeM])),cmap='gray',interpolation="none", extent = spatialExtent)
    #plt.show()
    
    # First Fourier Plane
    plt.figure(3)
    ax=plt.axes()
    ax.set_title("$First$ $Fourier$ $Plane,$ $z=f_{1}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x} [mm]$",labelpad=8)
    ax.set_ylabel("$N_{y} [mm]$")
    plt.imshow((np.abs(np.fft.ifftshift(SLMInput)[viewRangeN:viewRangeM,viewRangeN:viewRangeM])**2),cmap='gray',interpolation="none", extent = spatialExtent)
    cb=plt.colorbar()
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Helicoidal Filter Mask
    plt.figure(4)
    ax=plt.axes()
    ax.set_title("$SPP.$ $m=%d,$ $N=%d,$ $z=f_{1}$"%(Lvor,NG),fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x} [mm]$",labelpad=8)
    ax.set_ylabel("$N_{y} [mm]$")
    plt.imshow(np.fft.fftshift(np.angle(SLMfilterMask)),cmap='gray',interpolation="none", extent = spatialExtent)
    cb=plt.colorbar()
    cb.set_label(r"$Phase$ $value$",fontsize=12)
    #plt.show()
    
    # Signal with Filter at Modulating Plane
    fig=plt.figure(figsize=(14,8))
    
    ax0=plt.subplot(1,2,1)
    ax0.set_title("$Field$ $phase$ $after$ $SPP,$ $z=f_{1}$",fontsize=14,position=(0.5,1.0))
    ax0.set_xlabel("$N_{x} [mm]$",labelpad=8)
    ax0.set_ylabel("$N_{y} [mm]$") 
    map0=ax0.imshow(np.angle(np.fft.ifftshift(SLMPlane)[viewRangeN:viewRangeM,viewRangeN:viewRangeM]),cmap='gray',interpolation="none", extent = spatialExtent)
    cb=plt.colorbar(map0,orientation='horizontal')
    cb.set_label(r"$Phase$ $value$",fontsize=12)
    
    ax1=plt.subplot(1,2,2)
    ax1.set_title("$Intensity$ $field$ $after$ $SPP,$ $z=f_{1}$",fontsize=14,position=(0.5,1.0))
    ax1.set_xlabel("$N_{x} [mm]$",labelpad=8)
    ax1.set_ylabel("$N_{y} [mm]$") 
    map1=ax1.imshow((np.abs(np.fft.fftshift(SLMPlane)[viewRangeN:viewRangeM,viewRangeN:viewRangeM])**2),interpolation="none", extent = spatialExtent)
    cb=plt.colorbar(map1,orientation='horizontal')
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Lyot Plane - Filtered Signal - Inverse Fourier Field
    plt.figure(6)
    ax=plt.axes()
    ax.set_title("$Field$ $at$ $z=f_{1}+2*f_{2}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x} [mm]$",labelpad=8)
    ax.set_ylabel("$N_{y} [mm]$") 
    plt.imshow(np.abs(lyotPlane[viewRangeN:viewRangeM,viewRangeN:viewRangeM])**2,cmap='gray',interpolation="none", extent = spatialExtent)
    cb=plt.colorbar()
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Lyod's Aperture - Signal truncation
    plt.figure(7)
    ax = plt.axes()
    ax.set_title("$Field$ $*$ $App.Stop,$ $z=f_{1}+2*f_{2}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$ [mm]",labelpad=8)
    ax.set_ylabel("$N_{y}$ [mm]") 
    plt.imshow((np.abs(lyotAperturePlane[viewRangeN:viewRangeM,viewRangeN:viewRangeM])**2),cmap='gray',interpolation="none", extent = spatialExtent)
    cb=plt.colorbar(orientation='vertical')
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Focalized plane region - PSF
    
    plt.figure(8)
    ax = plt.axes()
    ax.set_title('Intensity - PSF Stellar Leakage at camera plane: 4F',fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$ [mm]",labelpad=8,size=20)
    ax.set_ylabel("$N_{y}$ [mm]",size=20) 
    plt.imshow((np.abs((np.fft.ifftshift(outputField))[viewRangeN:viewRangeM,viewRangeN:viewRangeM])**2),interpolation="none", extent = spatialExtent)#,cmap='gray')
    cb=plt.colorbar(orientation='vertical')
    cb.set_label(r"$Intensity$ $units$",fontsize=12)

    plt.show()