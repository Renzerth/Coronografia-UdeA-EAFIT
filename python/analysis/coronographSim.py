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

#%%--------------
#PROGRAM SETTINGS
#----------------
        
plotsEnabled = True
#%%---------------
#System Parameters
#-----------------
        
Lvor = 3 # Topological Charge
NG = 256 # Number of gray levels
spatialSampling = 60.1e-3 # SLM Pixel Pitch (mm)
apertureRadius = 2.0 # Telescope - Lyot plane (mm)
#%%-------------------
#Vortex Analyzer Tools
#---------------------

vortexTools = vortexProfiler(dx=spatialSampling,radius=apertureRadius)
#%%----------------------
#Compute Field Properties
#------------------------

aperture = np.fft.fftshift(vortexTools.placeAperture())
SLMfilterMask = np.fft.fftshift(vortexTools.SPP(Lvor,NG))
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
    
    scaleRange = 0.2
    pixelShift = 1-(vortexTools.spaceSamples % 2) # Center graph if even matrix size is used
    
    viewRangeN = int((1 - scaleRange)*vortexTools.halfSamples)+ pixelShift
    viewRangeM = int((1 + scaleRange)*vortexTools.halfSamples) + pixelShift
    
    # Truncated Input Plane at Stop
    plt.figure(2)
    ax=plt.axes()
    ax.set_title("$Telescope$ $Apperture,$ $z=0$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$")
    plt.imshow((abs(np.fft.ifftshift(aperture)[viewRangeN:viewRangeM,viewRangeN:viewRangeM])**2),cmap='gray')
    #plt.show()
    
    # First Fourier Plane
    plt.figure(3)
    ax=plt.axes()
    ax.set_title("$First$ $Fourier$ $Plane,$ $z=f_{2}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$")
    plt.imshow((abs(np.fft.ifftshift(SLMInput)[viewRangeN:viewRangeM,viewRangeN:viewRangeM])**2),cmap='gray')
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
    map0=ax0.imshow(np.angle(np.fft.ifftshift(SLMPlane)[viewRangeN:viewRangeM,viewRangeN:viewRangeM]),cmap='gray')
    cb=plt.colorbar(map0,orientation='horizontal')
    cb.set_label(r"$Phase$ $value$",fontsize=12)
    
    ax1=plt.subplot(1,2,2)
    ax1.set_title("$Intensity$ $field$ $after$ $SPP,$ $z=f_{1}$",fontsize=14,position=(0.5,1.0))
    ax1.set_xlabel("$N_{x}$",labelpad=8)
    ax1.set_ylabel("$N_{y}$") 
    map1=ax1.imshow((abs(np.fft.fftshift(SLMPlane)[viewRangeN:viewRangeM,viewRangeN:viewRangeM])**2))
    cb=plt.colorbar(map1,orientation='horizontal')
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Lyot Plane - Filtered Signal - Inverse Fourier Field
    plt.figure(6)
    ax=plt.axes()
    ax.set_title("$Field$ $at$ $z=f_{1}+2*f_{2}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$") 
    plt.imshow((abs(np.fft.ifftshift(lyotPlane)[viewRangeN:viewRangeM,viewRangeN:viewRangeM])**2),cmap='gray')
    cb=plt.colorbar()
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Lyod's Aperture - Signal truncation
    plt.figure(7)
    ax = plt.axes()
    ax.set_title("$Field$ $*$ $App.Stop,$ $z=f_{1}+2*f_{2}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$") 
    plt.imshow((abs(np.fft.ifftshift(lyotAperturePlane)[viewRangeN:viewRangeM,viewRangeN:viewRangeM])**2),cmap='gray')
    cb=plt.colorbar(orientation='vertical')
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Focalized plane region - PSF
    
    plt.figure(8)
    ax = plt.axes()
    ax.set_title("$Intensity$ $field$ $at$ $camera,$ $z=f_{1}+2*f_{2}+2f_{3}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$") 
    plt.imshow((abs(20*np.log10(np.fft.ifftshift(outputField))[viewRangeN:viewRangeM,viewRangeN:viewRangeM])**2))#,cmap='gray')
    cb=plt.colorbar(orientation='vertical')
    cb.set_label(r"$Intensity$ $units$",fontsize=12)

    plt.show()