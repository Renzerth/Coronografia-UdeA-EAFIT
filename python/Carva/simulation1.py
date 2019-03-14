"""
Author: Juan Pablo Carvajal.
Last update in 29-01-2019.
For GOF - UdeA.

Programm description: Basically, simulates the propagation of a Gaussian beam through a simple coronograph. 
Small variations have been done to simulate different situations. 
The current version was used to animate the evolution of the results while changing the topological charge (fractionary and integer).
"""




import numpy as np
import matplotlib.pyplot as plt
import scipy.interpolate as itp

# PROPAGATION FUNCS

# TF propagation func

def propTF(u1,L,lambd,z):

    M = len(u1)

    dx = L/M   # L is source and observation plane side length

    k  = 2*np.pi/lambd

    fx = np.arange(-1/(2*dx),1/(2*dx),1/L)

    FX,FY = np.meshgrid(fx,fx)

    H = np.exp(-1j*np.pi*lambd*z*(FX**2+FY**2))  # Transfer function
    H = np.fft.fftshift(H)

    # source fourier transform

    U1 = np.fft.fft2(np.fft.fftshift(u1))

    # observer fouriered field (convolution th.)

    U2 = H*U1

    # obs field

    u2 = np.fft.ifftshift(np.fft.ifft2(U2))

    return u2

# IR propagation func

def propIR(u1,L,lambd,z):

    M, N = len(u1),len(u1)

    dx = L/M   # L is source and observation plane side length

    k  = 2*np.pi/lambd

    x = np.arange(-L/2,L/2,dx)

    X,Y = np.meshgrid(x,x)

    h = 1/(1j*lambd*z)*np.exp(1j*k/(2*z)*(X**2+Y**2))  # Impulse response

    H = np.fft.fft2(np.fft.fftshift(h))*dx**2

    # source fourier transform

    U1 = np.fft.fft2(np.fft.fftshift(u1))

    # observer fouriered field (convolution th.)

    U2 = H*U1

    # obs field

    u2 = np.fft.ifftshift(np.fft.ifft2(U2))

    return u2

# OPTICAL ELEMENTS FUNCTIONS

# field through lens func

def lens(U, X, Y, f_lens, r_lens, wl):
    
    r = np.sqrt(X**2 + Y**2)
    LF = np.exp(-(1j*np.pi)/(wl*f_lens)*r**2)
    
    return U*LF
    
# field through cylindrical lense (x oriented) func

def clens(U, X, Y, f_lens, r_lens, wl): # wl: wavelenght
    
    LF = np.exp(-(1j*np.pi)/(wl*f_lens)*X**2)
    
# continuous spp pass-by

def sppc(U,X,Y,m):
    theta = np.arctan2(X,Y)
    vortex = np.exp(-1j*m*theta)
    return vortex*U

# discrete spp pass-by

def sppd(U, X, Y, m, l):
    theta = 2*np.pi/l*np.around(l*np.arctan2(X,Y)/(2*np.pi))
    vortex = np.exp(-1j*m*theta)
    
    return vortex*U

# SYSTEM PROPIERTIES

# entrance stop
r_e = 0.015 #meters

# lens 1
r_l1 = r_e
f_l1 = 0.3 #meters

# SPP
spp = []
for i in range(41):
    spp.append([float(i)/20. + 1.,256])

# Generating some field.

# U1 is the field that entries to the stop entry, for example, a gaussian.

# some grid values in a 1m square

L = 0.04 #
M = 4000 # for memory reason, use less than 4000. It will take some time for high numbers
dx = L/M
lambda_ = 0.633e-6 # wavelength (cm)

x, y = np.arange(-L/2,L/2,dx), np.arange(-L/2,L/2,dx)
X,Y = np.meshgrid(x,y)

# Gaussian func. for the grid

a = 500000 # some const. for the field

U1 = np.exp(-a*(X**2+Y**2))

# Now if we want a circular stop or radiud r = 0.5, is as simple as:

UC = np.sqrt(X**2+Y**2) < r_l1 # boolean-like array for simplifying the calculations

# Entry stop
U1 = U1*UC # U1 is the field at entry stop


# ENTIRE OPTICAL SYSTEM FUNCTION FOR AN INITIAL FIELD AND DETERMINED SPP

def os(U1, sppm):
    # lens1
    # lens(U, X, Y, f_lens, r_lens, wl)
    U1 = lens(U1, X, Y, f_l1, r_l1, lambda_)*UC # post lens1

    # prop1
    # propTF(u1,L,lambd,z)
    U1 = propTF(U1, L, lambda_, f_l1)*UC # pre spp

    # spp1
    # sppd(U, X, Y, m, l)
    U1 = sppd(U1,X,Y,sppm[0],sppm[1])*UC # spp
    
    # vortex lens
    U1 = clens(U1,X,Y,f_l1,r_l1, lambda_)*UC # cilyndrical lens

    # prop to lens 2
    # propTF(u1,L,lambd,z)
    U1 = propTF(U1, L, lambda_, f_l1)*UC # post spp

    # lens2
    # lens(U, X, Y, f_lens, r_lens, wl)
    U1 = lens(U1, X, Y, f_l1, r_l1, lambda_)*UC

    # prop to output stop
    U1 = propTF(U1,L,lambda_,f_l1)*UC
    
    return U1


for i in range(41):

    U2 = os(U1,spp[i])

    dis = '-varm-'+str(i)

    plt.imshow(np.abs(U2)**2)
    plt.title('m={:.2f}'.format(spp[i][0]))
    plt.colorbar()
    plt.savefig('mag'+dis+'.png')

    plt.clf()

    plt.imshow(np.angle(U2))
    plt.title('m={:.2f}'.format(spp[i][0]))
    plt.colorbar()    

    plt.savefig('phase'+dis+'.png')
    print(i)
    plt.clf()
