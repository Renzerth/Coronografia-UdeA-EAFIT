# -*- coding: utf-8 -*-
"""
Created on Tue Jul 28 17:52:06 2015

@author: elgaral

Last version: Edgar Rueda, october 2017

"""

import numpy as np
import cmath as cm
import scipy.misc

s=[]

"""
###### Discretization ####

"""
def Levels(Ima,L,tipo='medio'):
    maxval = np.max(Ima)
    Ima = L*(Ima/np.max(Ima))
    if tipo == 'piso':
        Ima = np.floor(Ima)
    elif tipo == 'techo':
        Ima = np.ceil(Ima)
    elif tipo == 'zp':
        Ima = 2*(Ima/np.max(Ima))
        Ima = abs(np.floor(Ima) - 1)
        return maxval*Ima
    else: # 'medio'
        Ima = np.floor(Ima) + 0.5

    return maxval*(Ima/L)



"""
###### Astigmatism #####

N: matrix size
dx: sample size
wl:  wavelength
f:   focal distance

Lens: Lens complex amplitude
"""
def Astig(N,dx):
    x = np.arange(-N*dx/2,N*dx/2,dx)  # -N*dx/2 <= x <= N*dx/2 - dx
    X , Y = np.meshgrid(x,x)   #
    Fase = np.exp(1j*1.2*(X**2 + Y**2)*np.sin(2.0*np.arctan2(Y,X)))
    return Fase



"""
###### Circular pupil #######
Creates a circular pupil inside a matrix

a,b: index of center of matrix (in real matrix position)
N: base array size
radius: pupil radius in pixels

array: final matrix
"""
def cmask(a,b,radius,N):
  y,x = np.ogrid[-a:N-a,-b:N-b]
  mask = x*x + y*y <= radius*radius
  array = np.zeros((N, N))
  array[mask] = 1.

  return array

"""
###### Embed a matrix into other #######
Embeded a matrix into other matrix of bigger size. The sizes must be multiples

base = base matrix
center = matrix to be embed

final: final matrix
"""
def matEmbed(base,center):
    N = int(np.sqrt(base.size))
    NN = int(np.sqrt(center.size))
    base[int((N-NN)/2):int((N+NN)/2),int((N-NN)/2):int((N+NN)/2)] = center
    return base




"""
###### RGB t Gray ####
Converts a RGB image (array) into a grayscale image, using the luminosity formula

A : original matrix, RGB


imagen: final matrix, gray
"""
def rgb2gray(A):
    imagen = (0.21*A[:,:,0] + 0.72*A[:,:,1] + 0.07*A[:,:,2])
    return imagen

"""
###### zoom ####
Makes a zoom of the central part of an even matrix

matriz : original matrix
N: original size
NN: zoom size

zoom: final matrix
"""
def zoomC(N,NN,matriz):
    zoom = matriz[int((N-NN)/2):int((N-NN)/2)+NN,int((N-NN)/2):int((N-NN)/2)+NN]
    return zoom

"""
###### plane wave #####

N: matrix size
dx: sample size
wl:  wavelength
A0: plane wave amplitude magnitude
alfax:   director cosine angle with respect to x-axis in rads
alfay: director cosine agle with respect to y-axis in rads

Planew: Plane wave complex amplitude
"""
def planew(N,dx,alfax,alfay,wl,A0):
    k = 2*np.pi/wl  # wave number

    x = np.arange(-N*dx/2,N*dx/2,dx)  # -N*dx/2 <= x <= N*dx/2 - dx
    X , Y = np.meshgrid(x,x)   # coord plane
    Planew = np.exp(1j*k*(np.cos(alfax)*X + np.cos(alfay)*Y))
    return Planew

"""
###### Gaussian beam (Saleh, 3. 1-7) #####

N: matrix size
dx: sample size
W0:  beam waist radius
A0:  Beam maximum amplitud at waist
wl:  wavelength
z:   distance from the waist
info: 1 if you want the main beam information printed

G: Gaussian beam complex amplitude
"""
def Gbeam(N,dx,W0,A0,wl,z,info):
    if z == 0:
        z = 1.0e-20

    z0 = np.pi*W0**2/wl  # Rayleigh range
    W = W0*np.sqrt(1 + (z/z0)**2)  # Beam width at z
    k = 2*np.pi/wl # wave number
    R = z*(1 + (z0/z)**2) # beam radius of curvature
    Xi = np.arctan(z/z0) # Guoy phase

    x = np.arange(-N*dx/2,N*dx/2,dx)  # -N*dx/2 <= x <= N*dx/2 - dx
    X , Y = np.meshgrid(x,x)   #
    Amp = np.exp(-(X**2 + Y**2)/W**2)
    Phase = np.exp(-1j*k*z - 1j*k*(X**2 + Y**2)/(2*R) + 1j*Xi)

    G = A0*(W0/W)*Amp*Phase
    return G

"""
###### Thin Lens (Goodman, 5. 5-10) #####

N: matrix size
dx: sample size
wl:  wavelength
f:   focal distance

Lens: Lens complex amplitude
"""
def Lente(N,dx,f,wl):
    k = 2*np.pi/wl

    x = np.arange(-N*dx/2,N*dx/2,dx)  # -N*dx/2 <= x <= N*dx/2 - dx
    X , Y = np.meshgrid(x,x)   #
    Lens = np.exp(-1j*k/(2*f)*(X**2 + Y**2))
    return Lens

"""
###### Phase of a complex Image #####

N: matrix size
Ima:  Complex Image

Fase: Phase of the complex image
"""
def FaseIma(N,Ima):
    Fase = np.zeros((N,N))
    for cont1 in np.arange(0,N):
        for cont2 in np.arange(0,N):
            Fase[cont1,cont2] = cm.phase(Ima[cont1,cont2])
    return Fase


"""
###### SPP  #####

N: matrix size
dx: sample size
ch: SPP charge

SPP: Spiral phase plate complex value
"""
def SPP(N,dx,ch):
    x = np.arange(-N*dx/2,N*dx/2,dx)  # -N*dx/2 <= x <= N*dx/2 - dx
    X , Y = np.meshgrid(x,x)   #

    angle = np.arctan2(Y,X)
    SPP = np.exp(1j*ch*angle)
    return SPP

"""
###### Scale Image #####

N: matrix size
Esc: scaling factor
Type: type of image, real amplitude (amp) or phase (fase)
G: Image of real values

G2: scaled image
"""
def scalaIma(N,Esc,Type,G):
    if Esc != 1:
        if Type == 'amp':
            maxi = np.amax(G)
            G = 255*G/maxi
            G = scipy.misc.imresize((G),Esc)
            G = G*maxi/255.
            print(maxi)
        elif Type == 'fase':
            G = 255*(G + np.pi)/(2*np.pi)
            G = scipy.misc.imresize(G,Esc)
            G = 2*np.pi*G/255. - np.pi


        NN = np.round(N*Esc)
        if NN % 2 != 0: #odd
            NN = NN  -1

        if NN > N:
            G2 = G[(NN-N)/2:(NN-N)/2+N,(NN-N)/2:(NN-N)/2+N]

        if NN < N:
            G2 = np.zeros((N,N))
            G2[(N-NN)/2:(N-NN)/2+NN,(N-NN)/2:(N-NN)/2+NN] = G
    else:
        G2 = G

    return G2



def clear_all():
    """Clears all the variables from the workspace of the spyder application."""
    gl = globals().copy()
    for var in gl:
        if var[0] == '_': continue
        if 'func' in str(globals()[var]): continue
        if 'module' in str(globals()[var]): continue

        del globals()[var]
if __name__ == "__main__":
    clear_all()
