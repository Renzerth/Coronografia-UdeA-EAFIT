#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Author: Carva.
Last update in 13-03-2019.
Coronografía UdeA-EAFIT.

Propagation funtions in functional form.


"""


import numpy as np


# TF propagation func

def propTF(u1,dx,wl,z):
    """
    Propagation - Transfer Function (TF) approach asumes
    same x and y side lengths and uniform sampling
    Inputs:
     -u1: source plane field
     -dx: source and observation plane side length delta
     -wl: wavelength
     -z: propagation distance (optical axis)
    Output:
     -us: observation plane field
    
    Chapter 5 - Computational Fourier Optics
    
    Explanation:
          g0 = fftshift(g); % Shift
          G0 = fft2(g0)*dA; % 2D fft and dxdy scaling (dx = dy)
          G = fftshift(G0); % Center
    """

    M = len(u1)

    L = dx/M   # L is source and observation plane side length

    #k  = 2*np.pi/wl

    fx = np.arange(-1/(2*dx),1/(2*dx),1/L)

    FX,FY = np.meshgrid(fx,fx)

    H = np.exp(-1j*np.pi*wl*z*(FX**2+FY**2))  # Transfer function
    H = np.fft.fftshift(H)

    # source fourier transform

    U1 = np.fft.fft2(np.fft.fftshift(u1))

    # observer fouriered field (convolution th.)

    U2 = H*U1 # at this point both H and U1 are shifted, to reduce operations and optimize

    # obs field

    u2 = np.fft.ifftshift(np.fft.ifft2(U2))

    return u2

# IR propagation func

def propIR(u1,dx,wl,z):

    """
    Propagation - Impulse Response (IR) approach asumes
    same x and y side lengths and uniform sampling
    Inputs:
     -u1: source plane field
     -dx: source and observation plane side length delta
     -wl: wavelength
     -z: propagation distance (optical axis)
    Output:
     -us: observation plane field
    
    Chapter 5 - Computational Fourier Optics
    
    Explanation:
          g0 = fftshift(g); % Shift
          G0 = fft2(g0)*dA; % 2D fft and dxdy scaling (dx = dy)
          G = fftshift(G0); % Center
    """

    M = len(u1)

    L = dx/M   # L is source and observation plane side length

    k  = 2*np.pi/wl

    x = np.arange(-L/2,L/2,dx)

    X,Y = np.meshgrid(x,x)

    h = 1/(1j*wl*z)*np.exp(1j*k/(2*z)*(X**2+Y**2))  # Impulse response

    H = np.fft.fft2(np.fft.fftshift(h))*dx**2

    # source fourier transform

    U1 = np.fft.fft2(np.fft.fftshift(u1))

    # observer fouriered field (convolution th.)

    U2 = H*U1

    # obs field

    u2 = np.fft.ifftshift(np.fft.ifft2(U2))

    return u2

# FF propagation func.

def propFF(u1,dx1,wl,z):

    """
    propagation - FraunhoFer (FF) pattern assumes uniform sampling
    Inputs:
     -u1: source plane field
     -dx1: source plane side length delta
     -wl: wavelength
     -z: propagation distance
    Outputs:
     -L2: observation plane side length
     -u2: observation plane field
    
    Chapter 5 - Computational Fourier Optics
        --Based on Saples matlab work."""
    
    M = len(u1)
    
    L1 = dx1/M
    
    k = 2*np.pi/wl
    
    L2 = wl*z/dx1
    dx2 = wl*z/L1
    x2 = np.arange(-L2/2,L2/2,dx2)
    
    X2,Y2 = np.meshgrid(x2,x2)
    
    c = 1/(1j*wl*z)*np.exp(1j*k/(2*z)*(X2**2+Y2**2))
    u2 = c*np.fft.ifftshift(np.fft.fft2(np.fft.fftshift(u1))) * dx1**2
    
    return u2, L2
    
# AS propagation (angular spectrum) func.

def propAS(u1, dx, wl, z):

    """
    This function makes propagation of an input field based on the angular
    spectrum of plane waves presented in chapter 3.10 of Goodman.
    
    Assumes uniform sampling and squared planes
    
    Input :
     
    U1 : Input field
    dx : source and observation planes side lenght delta
    wl : Ilumination wavelength
    z : Propagation distance
    
    Output :
     
    U2 : Output field
    
        --Based and modified from Saples matlab code.
    """
    
    M = len(u1)
    
    L = dx/M
    
    #x = np.arange(-L/2,L/2,dx)
    fx = np.arange(-1/(2*dx),1/(2*dx),1/L)
    
    fX, fY = np.meshgrid(fx,fx)
    
    #fR = np.sqrt(fX**2+fY**2)
    
    #mask = fR < 1/wl
    
    #k = 2*pi/wl
    
    A1 = np.fft.fftshift(np.fft.fft2(np.fft.fftshift(u1)))
    
    #A2 = A1*mask*np.exp(1j*k*z*np.sqrt(1-(wl*fX)**2-(wl*fY)**2))
    A2 = A1*np.exp(1j*k*z*np.sqrt(1-(wl*fX)**2-(wl*fY)**2))
    
    U2 = np.fft.fftshift(np.fft.ifft2(np.fft.ifftshift(A2)))
    
    return U2
    
    
# Translate to py3.6 the Sypek propagator in the future if necesary.
    
"""
##### Sypek propagator (OC, v.116, 1995) #####
u1:  input complex field
dx:  sample size
wl:  wavelength
z:   propagation distance

Ima3: exit complex field

################################################

def propSypek(u1,dx,wl,z):
    M = u1.shape[0]   # input field array dimension (square)
    L = M*dx # working side length
    
    zcr = L*dx/wl #propagation distance in mm, critical sampling
    steps = math.floor(z/zcr) # number of propagations to keep the critical sampling
    print 'Propagation distance z = %.2f mm' % z
    print 'Critical propagation distance z = %.2f mm' % zcr
    print 'Number of propagaton steps = %d' % steps
    
    if steps != 0:
        extra = z - zcr*steps # In the case that the total distance is not a multiple of the critical sampling distance
        print 'Extra propagation distance = %.2f mm' % extra
            
    N=2*M # array size for the fft2
    Ima3 = u1 # input field
    
    if steps == 0: # in the case that the propagation distance is smaller than the critical sampling distance
        Ima1 = np.zeros((N,N),'complex') # woking array
        Ima1[N/4:N/4+M+0,N/4:N/4+M+0] = Ima3 
        Ima2 = propTF(Ima1,dx,wl,z) # propagation
        Ima3 = Ima2[N/4:N/4+M+0,N/4:N/4+M+0] # exit field
        
    else: 
        for cont1 in np.arange(1,steps + 1): # propagation in steps PREGUNTA 0 o 1
            Ima1 = np.zeros((N,N),'complex') # OJO tiene que ser una matriz cumpleja
            Ima1[N/4:N/4+M+0,N/4:N/4+M+0] = Ima3
            Ima2 = propTF(Ima1,dx,wl,zcr)
            Ima3 = Ima2[N/4:N/4+M+0,N/4:N/4+M+0]
    
            
        if extra != 0.: # propagation of remaining distance
            Ima1 = np.zeros((N,N),'complex') # OJO tiene que ser una matriz cumpleja
            Ima1[N/4:N/4+M+0,N/4:N/4+M+0] = Ima3
            Ima2 = propTF(Ima1,dx,wl,extra)
            Ima3 = Ima2[N/4:N/4+M+0,N/4:N/4+M+0]
   
    #plt.imshow(abs(Ima3), cmap='gray',vmin=0.,vmax=1.,interpolation='nearest')
    #plt.colorbar(im, orientation='horizontal')
    #plt.axis('off')
    #plt.title('Exit field magnitude')
    #plt.savefig("test.png",bbox_inches='tight')
    #plt.show()
    
    return Ima3
    
    
    """
