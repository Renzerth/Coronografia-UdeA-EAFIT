# -*- coding: utf-8 -*-
"""
Created on Tue Jun  2 19:43:24 2015

@author: Edgar Rueda

PROPAGATORS
"""

# libraries
import numpy as np
import cmath as cm
import math
import pylab as plt

s=[]

"""
###### Transfer Function propagator (Diels) #####


u1:  input complex field
dx:  sample size
wl:  wavelength
z:   propagation distance

u2: exit complex field
"""

def propTF(u1,dx,wl,z):
    
    M = u1.shape[0]   # input field array dimension (square)
    #dx = L/M   # sample interval
    L = M*dx
    #k = 2*np.pi/wl   # wavenumber
    
    fx = np.arange(-1/(2*dx),1/(2*dx),1/L)  # -1/(2*dx) <= fx <= 1/(2*dx) - 1/L
    FX , FY = np.meshgrid(fx,fx)   #
    
    H = np.zeros((M,M))+1j*np.zeros((M,M))
    for cont1 in np.arange(0,M):
        for cont2 in np.arange(0,M):
            val = cm.exp(-1j*np.pi*wl*z*(FX[cont1,cont2]**2+FY[cont1,cont2]**2))
            #print val
            H[cont1,cont2] = val   # Transfer fucntion
            
            
    H = np.fft.fftshift(H)   # shift transfer function
    U1 = np.fft.fft2(np.fft.fftshift(u1))
    U2 = H*U1
    u2 = np.fft.fftshift(np.fft.ifft2(U2))
    
    return u2
    
"""
##### Sypek propagator (OC, v.116, 1995)####
u1:  input complex field
dx:  sample size
wl:  wavelength
z:   propagation distance

Ima3: exit complex field
"""
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
