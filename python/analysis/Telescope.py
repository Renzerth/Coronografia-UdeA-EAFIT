#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on 28-07-2016
by AFIC

Last modification: 30-09-2018

This code emulates the propagation of an incoming plane wave 
throughout a Keplerian telescope. A generalized inclination might be 
included on the plane wave in order to simulate astronomical objects 
shifted with respect to the main axis of the optical system.

Modded by Grupo de Optica Aplicada - EAFIT
Juan José Cadavid Muñoz - 23-01-2019
"""

import sys, os; 
import numpy as np
sys.path.insert(0, os.path.abspath('..'))

import tools.funcionesMOD as fM
import tools.coherImagLIB as CIL
import tools.vortexGEN as vGEN #Many levels
import matplotlib.pyplot as plt
#from matplotlib import gridspec
#%%-------------------------
#PROGRAM SETTINGS
#-------------------------
plotsEnabled = True

#%%-------------------------
#OPTICAL PARAMETERS (in mm)
#-------------------------
f1 = 200. # Objective lens focal length
f2 = 200. # Ocular lens focal length
R2 = 8. # Apperture and Lyot stop radius
f3 = 300. 
wl = 532.e-6 # wavelength in mm  (green)

#---------------------------
#SPP: CHARGE AND GRAY LEVELS
#---------------------------
Lvor = 2 # Topological Charge
NG = 256 # Number of gray levels

#---------
#SAMPLING
#---------
M = 2**12 #2**12=4096 # Number of samples
dx = 26.e-3 # Modulator's pixel size in mm
L = M*dx #100. # Plane size in mm

#------------
#Data Slicing
#------------
hM = int(M/2.0) # Half M

#------------------
#OPTICAL PARAMETERS
#------------------
Factor = 1/L/dx
k = 2.*np.pi/wl

Ray = 1.22*wl/(2*R2)
RayL = Ray*f2 #Approx. size of Airy disk at plane f2: Tan(Ray)=RayL/f2-->Ray=RayL/f2
RealS = 6/L*wl*f2 #At radius of Airy disk there are approx. 6 samples,df=1/L,scale=wl*f2  
                  # RayL is approx RealS
                  #  Computer Airy size: AiryC = 6samples*dx
                  #   RayL = 6*dx*Factor*wl*f2
                  #    thus, AiryC = RayL/(Factor*wl*f2)

AiryC = RayL/(Factor*wl*f2)

#Neptune: 4.3e9km from Earth (mV=8) 
# Triton: 3.6e5km from Neptune (mV=13) if Delta_mV=5 --> I1=2.5**5*I2=100*I2
#  Nep-Tri: alpha from Earth is approx. 2*Ray <--> Tan(alpha)=3.6e5/4.3e9 
#   Nereid: 5.5e6km from Neptune (mV=19 --> 1e5 times less brilliant than Nep.)
#    Nep-Ner: alpha from Earth is approx. 22*Ray 
#     Pluto: (mV=15) Delta_mV=7 --> I1=2.5**7*I3=610*I3

alpha = 0*2*Ray #1st angular separation in rad
phi = 0*4*Ray #2nd angular separation in rad

#%%------------
#PLANET OFFSET
#------------

Roffx = np.tan(alpha)*f2 #Real optical x-axis offset in mm
Roffy = np.tan(phi)*f2   #y-axis...
Roff = np.sqrt(Roffx**2+Roffy**2)

Coffx = Roffx/(Factor*wl*f2) #Computational x-axis offset in mm 
Coffy = Roffy/(Factor*wl*f2) #y-axis...
Coff = np.sqrt(Coffx**2+Coffy**2) #=np.sqrt((OffsetSamples*dx)**2+(OffsetSamples*dy)**2)

#If alpha=0.00012 --> the real offset is approx. alpha*f2=0.036, the same than
 #22*1/L*wl*f2=0.035 --> 22 offset samples in frequency(1/L) with the scale wl*f2,
  #the same than 22*dx*Factor*wl*f2, where Factor is the relation between df and dx
   #Factor=(1/L)/dx. The offset given by the computer is Coff=22*dx=Roff/(Factor*wl*f2)   

#To reach AiryC we need an AiryAngle such that AiryC=Coffx, Roffx=RayL !! we already have
# the condition to reach AiryC: we need to reach in the reality RayL, which is reached
#  with Ray, of this way due to the explanation of the last paragraph we get AiryC:
#   with Ray we get RayL=Roffx in the real life but in the computer we get AiryC=Coffx. 

#--------------
#APPERTURE STOP
#--------------
x=np.arange(-L/2,L/2,dx)
y=np.arange(-L/2,L/2,dx)

AS = np.zeros((M,M))
lim1 = hM - (int(R2/dx) + 10) #Number of pixels needed to reach R2 plus a 
lim2 = hM + (int(R2/dx) + 10) # little extra-range, centered in M/2 (This process is made to avoid big calculations with the whole MxM matrix)
AS[lim1:lim2,lim1:lim2] = fM.circ(x[lim1],x[lim2],dx,2*R2,0)    

                                     #2*R2: Function expects the diameter
                                      #does not return values of x and 
                                       #gives the circular shape of the App.Stop

#%%--------------
#INITIAL FIELD
#--------------
z = 0
X,Y = np.meshgrid(x,y)
E0 = np.ones((M,M),'complex')
#Let's include the phase term due to the inclination of the input plane:
E0 = E0*np.exp(1j*k*(z*np.cos(alpha)*np.cos(phi)+X*np.sin(alpha)+
                     Y*np.cos(alpha)*np.sin(phi)))

#If one has a list rather than an array, it's necessary to work with the notation L[][] to 
# call an element from it, this is the case of X (see above), which is a natural list

Fase = CIL.phase(E0)

#%%-------------------------------
#FIELD AFTER THE FIRST APPERTURE
#-------------------------------
E1 = 0 # No binary body is available
E0 = (E0+E1)*AS # Truncated Input field

dx2 = dx**2
print("==================================================")
print('Initial Field Power', dx2*np.sum(abs(E0)**2))

#--------------------------
#FIRST FOURIER PLANE (z=f1)
#--------------------------
#At the focal plane we have the Fourier transform of the initial Field times the pupil function

fx2 = np.arange(-1/(2.*dx),1/(2.*dx),1./L)
fy2 = fx2
 
Fx2,Fy2 = np.meshgrid(fx2,fy2)

Ef2 = np.fft.fftshift(E0)
Ef2 = np.fft.fft2(Ef2)
Ef2 = np.fft.ifftshift(Ef2)
Ef2 = Ef2*dx*dx
# (HAS TO BE SCALED by an additional phase factor when the field is 
#  just against the lens (d = 0) and not at d = -f1; see GOODMAN_ed2 eq 5.16)
# Ef2 = Ef2* np.exp(1j*k*f1*wl*wl*(Fx2*Fx2+Fy2*Fy2))#/(1j*wl*f1)

print("Incoming max Intensity = ",(abs(Ef2)**2).max()) #Pot/m**2. In a plane wave I=c*epsilon*E**2-->
                                                      # See Hector_Alzate Pag.128 Fisica de las ondas
print("Field at Stop - Power:", L**-2*np.sum(abs(Ef2)**2))         #Power=I*m**2 

#%%----------------------
#INTRODUCING SPP (z=f1)
#----------------------
#Now, let's model the SPP using its characteristic transmitance times a circ function  
C = np.zeros((M,M),'complex')
C = vGEN.SPP(Lvor,NG,L,dx)

FaseC = CIL.phase(C)

#%%---------------------
#FIELD PASSES BY C
#---------------------
Ev = Ef2*C
FaseEv=CIL.phase(Ev)

print("Field at SLM Plane - Power:", L**-2*np.sum(abs(Ev)**2))

#%%------------------
#FIELD AT z=f1+2*f2
#------------------
EvP = np.fft.fftshift(Ev)
EvP = np.fft.fft2(EvP)
EvP = np.fft.ifftshift(EvP)
EvP = EvP*1./L**2#/M#*dx*dx
#EvP = EvP #/(1j*wl*f2)

print("Lyot's Plane Distribution - Power:", dx2*np.sum(abs(EvP)**2))#*1./L**2)

#----------
#FILTERING
#----------
Rf = R2
Aff = np.zeros((M,M))
lim1f = int(M/2 - (Rf/dx + 10)) #Number of pixels needed to reach Rf plus a                                                                                                               
lim2f = int(M/2 + (Rf/dx + 10)) # little extra-range, centered in M/2      
                                                                                                               
Aff[lim1f:lim2f,lim1f:lim2f] = fM.circ(x[lim1f],x[lim2f],dx,2*Rf,0)
Ef = EvP*Aff

print("Lyot's Aperture Truncated Field - Power", dx2*np.sum(abs(Ef)**2))
print("==================================================")

#-----------
#FINAL FIELD
#-----------
Eff  = np.fft.fftshift(Ef)
Eff = np.fft.fft2(Eff)
Eff = np.fft.ifftshift(Eff)
Eff = Eff*dx*dx
Eff = Eff #*np.exp(1j*k*f3*wl*wl*(Fx2*Fx2+Fy2*Fy2))#/(1j*wl*f3)

#%%---
#PLOTS
#-----
if plotsEnabled:
    plt.close()
    
    # Visualization Ranges
    viewRangeNA = hM - 800
    viewRangeMA = hM + 800
    
    viewRangeNB = hM - 148
    viewRangeMB = hM + 148
    
    viewRangeNC = hM - 120
    viewRangeMC = hM + 120
    
    # Input Plane
    plt.figure(1)
    ax=plt.axes()
    ax.set_title("$Initial$ $Field$ $phase$")
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$")
    plt.imshow(Fase)
    cb=plt.colorbar()
    cb.set_label("$Phase$ $Value$")
    #plt.show()
    
    # Truncated Input Plane at Stop
    plt.figure(2)
    ax=plt.axes()
    ax.set_title("$Telescope$ $Apperture,$ $z=0$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$")
    plt.imshow((abs(E0)**2)[viewRangeNA:viewRangeMA,viewRangeNA:viewRangeMA],cmap='gray')
    #plt.show()
    
    # First Fourier Plane
    plt.figure(3)
    ax=plt.axes()
    ax.set_title("$First$ $Fourier$ $Plane,$ $z=f_{2}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$")
    plt.imshow((abs(Ef2)**2)[viewRangeNB:viewRangeMB,viewRangeNB:viewRangeMB],cmap='gray')
    cb=plt.colorbar()
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Helicoidal Filter Mask
    plt.figure(4)
    ax=plt.axes()
    ax.set_title("$SPP.$ $m=%d,$ $N=%d,$ $z=f_{1}$"%(Lvor,NG),fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$")
    plt.imshow(FaseC,cmap='gray')
    cb=plt.colorbar()
    cb.set_label(r"$Phase$ $value$",fontsize=12)
    #plt.show()
    
    # Signal with Filter at Modulating Plane
    fig=plt.figure(figsize=(14,8))
    
    ax0=plt.subplot(1,2,1)
    ax0.set_title("$Field$ $phase$ $after$ $SPP,$ $z=f_{1}$",fontsize=14,position=(0.5,1.0))
    ax0.set_xlabel("$N_{x}$",labelpad=8)
    ax0.set_ylabel("$N_{y}$") 
    map0=ax0.imshow(FaseEv[viewRangeNC:viewRangeMC,viewRangeNC:viewRangeMC],cmap='gray')
    cb=plt.colorbar(map0,orientation='horizontal')
    cb.set_label(r"$Phase$ $value$",fontsize=12)
    
    ax1=plt.subplot(1,2,2)
    ax1.set_title("$Intensity$ $field$ $after$ $SPP,$ $z=f_{1}$",fontsize=14,position=(0.5,1.0))
    ax1.set_xlabel("$N_{x}$",labelpad=8)
    ax1.set_ylabel("$N_{y}$") 
    map1=ax1.imshow((abs(Ev)**2)[viewRangeNC:viewRangeMC,viewRangeNC:viewRangeMC])
    cb=plt.colorbar(map1,orientation='horizontal')
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Lyot Plane - Filtered Signal - Inverse Fourier Field
    plt.figure(6)
    ax=plt.axes()
    ax.set_title("$Field$ $at$ $z=f_{1}+2*f_{2}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$") 
    plt.imshow((abs(EvP)**2)[viewRangeNA:viewRangeMA,viewRangeNA:viewRangeMA],cmap='gray')
    cb=plt.colorbar()
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Lyod's Aperture - Signal truncation
    plt.figure(7)
    ax = plt.axes()
    ax.set_title("$Field$ $*$ $App.Stop,$ $z=f_{1}+2*f_{2}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$") 
    plt.imshow((abs(Ef)**2)[viewRangeNA:viewRangeMA,viewRangeNA:viewRangeMA],cmap='gray')
    cb=plt.colorbar(orientation='vertical')
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    #plt.show()
    
    # Focalized plane region - PSF
    
    plt.figure(8)
    ax = plt.axes()
    ax.set_title("$Intensity$ $field$ $at$ $camera,$ $z=f_{1}+2*f_{2}+2f_{3}$",fontsize=14,position=(0.5,1.0))
    ax.set_xlabel("$N_{x}$",labelpad=8)
    ax.set_ylabel("$N_{y}$") 
    plt.imshow((abs(Eff)**2)[viewRangeNC:viewRangeMC,viewRangeNC:viewRangeMC])#,cmap='gray')
    cb=plt.colorbar(orientation='vertical')
    cb.set_label(r"$Intensity$ $units$",fontsize=12)
    
    plt.show()

print("==================================================")
print("Final field power:", L**-2*np.sum(abs(Eff)**2))
print("Final max Intensity = ",(abs(Eff)**2).max())
print("Topological charge:", Lvor, "- Gray levels:", NG)
print("Screen size (mm), Num. of samples, Step size (mm):",L,M,dx)
print("==================================================")
#plt.close()

#--------
#--------

#f = open("writing.txt",'w')
#for i in np.arange(0,len(x)): f.write("%f\n"%x[i])