#!/usr/bin/env python3
# -*- coding: utf-8 -*-
#Utility: coherent imaging library.
# See an example at ...AFIC/Fourier_comp/PROPYTHON

from tools import funcionesMOD as fM
import numpy as np
#import cmath as cm

def lens(p,f):
    q = 1./(1./f-1./p)
    Mag = q/p
    return q,Mag

def field2D(Ex,Ey):
    Ex,Ey = np.meshgrid(Ex,Ey)
    E = Ex*Ey
    return E

def phase(E):
    #F = np.zeros(E.shape)
    F = np.angle(E)
#    for i in np.arange(0,E.shape[0]):
#        for j in np.arange(0,E.shape[1]):
#            F[i,j] = cm.phase(E[i,j])

    return F
    
def phaseFT(E):  #Returns the phase of the field after Fourier transform
    Eshift = np.fft.fftshift(E)
    E = np.fft.fft2(Eshift)
    E = np.fft.ifftshift(E)
#    F = np.zeros(E.shape)
#    for i in np.arange(0,E.shape[0]):
#        for j in np.arange(0,E.shape[1]):
#            F[i,j] = cm.phase(E[i,j])
    
    F = np.angle(E)

    return F
    
def imaging(E,P):
    Eshift = np.fft.fftshift(E)
    E = np.fft.fft2(Eshift)
    E = np.fft.ifftshift(E)
    E = P*E
    E = np.fft.fftshift(E)
    E = np.fft.ifft2(E) 
    ET = np.fft.ifftshift(E)
    return ET


    
