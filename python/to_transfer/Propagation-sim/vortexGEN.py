"""
Created on 29-07-2015
by AFIC

Last modification: 05-11-2018
"""
import numpy as np

def cart2pol(x, y):
    rho = np.sqrt(x**2 + y**2)
    phi = np.arctan2(y, x)
    return(rho, phi)

#SLModulator GOF, pixel size: 26 micrometers, 26e-3 mm
def SPP(Lvor,NG,L,dx=26.e-3,B=False):
    
    M = L/dx
    N = M
    x = np.arange(-L/2, L/2, dx)
    xx, yy = np.meshgrid(x, x)
    _,phi = cart2pol(xx, yy)

    phi1 = phi + np.pi #Matrix with entries within [0:2pi-step]
    phaseVor = np.mod(Lvor * phi1, 2*np.pi) #256-levels discretization
    
    phi2 = np.floor(phaseVor/(2*np.pi/NG)) #Matrix with whole-numbers between 0 and NG-1 
    phi3 = phi2/NG #phi3 is phi2 but normalized
    phi4 = np.round(256 * phi3) #256 fractions in each region
    
    return np.exp(1j*(2*np.pi*phi3 - np.pi))
