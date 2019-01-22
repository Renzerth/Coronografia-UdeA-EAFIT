# Prueba #
##########

from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
from scipy.special import gamma, j1, jv
from scipy.misc import factorial
from scipy.integrate import quad, simps

# parametros
N = 2
L = 1
jobs = 0
mobs = L*(1 + jobs*N)


f1 = 1.0
f2 = 0.2#0.1
fFR = 1.6
wl = 532e-9
z0 = f2 - mobs*f2**2 / (L * fFR) 


a = 5.e-2
b = 1000*a
w0 = a/10. 

k = 2*np.pi / wl
j = np.array([-1,0]) #np.array([-1,0])
m = L*(1 + j*N)
print (m)
#rmax = 2*10.e-6 #Gaussian case
rmax = a*f2/f1

NN = 63 #64/4
x = np.linspace(-2*rmax/1,2*rmax/1,NN)
#x = np.linspace(-5,5,NN)
X,Y = np.meshgrid(x,x)
R = np.sqrt(X**2+Y**2)
T = np.arctan2(Y,X)

def integrand_bessel(rho,m,r):
    alpha = 0.5*1j*k * (m/(L*fFR) + z0/f2**2 - 1/f1 - 1/f2)
    return j1(k*a/f1 * rho) * np.exp(-alpha * rho**2) * jv(m, k*r/f2 * rho)

def integrand_gauss(rho,m,r):
    #z0 = f2 - m*f2**2 / (L * fFR) 
    alpha = 0.5*1j*k * (m/(L*fFR) + z0/f2**2 - 1/f2) + 1/w0**2
    return 1j * np.exp(-alpha * rho**2) * jv(m, k*r/f2 * rho) * rho
    

def do_image(m, gauss = True):
    #rho = np.linspace(0,10*a,1000000)
    rho = np.linspace(0,10*a,1000000)
    R_unique = list(set(R.flatten()))
    u_m = np.zeros((NN,NN)) + 1j*np.zeros((NN,NN))

    arg = np.pi * m / (L*N)
    
    if gauss:
        integrand = integrand_gauss
        K1 = (np.exp(-1j*arg) * np.sinc(arg / np.pi) #The sinc in np is defined as sinc(x) = sin(pi*x) / pi*x 
              * (k/f2) * (1j**(3*m+1))
              * np.exp(1j*k*(f2+z0))
              * np.exp(1j*m*T))
    else:
        integrand = integrand_bessel
        K1 = (np.exp(-1j*arg) * np.sinc(arg / np.pi)
              * (k*a/f2) * (1j**(3*abs(m)-2))
              * np.exp(1j*k*(f2+z0))
              * np.exp(1j*m*T))

    id = 0
    
    
    for r in R_unique:
        I2 = np.trapz(integrand(rho, m, r), x = rho)
        indices = np.where(R == r)
        u_m[indices] = I2
        print (id)
        id+=1

    return K1 * u_m

U_RT = 0 #U(r,theta)
for i in m:
    U_RT += do_image(i, gauss = False)

INT = abs(U_RT)**2
INT_max = np.max(INT)
plt.imshow(INT / INT_max) #, norm = colors.LogNorm())
plt.colorbar()
plt.savefig("int_bessel._N%d_L%d_mobs%d.png"%(N,L,mobs))
#plt.savefig("int_gauss_N%d_L%d_mobs%d.png"%(N,L,mobs))
plt.show()

#integral from 0 to inf of (J1(2pi*a*rho / lamb*f1) * exp(-alpha*rho**2) * Jm(k*rho*r/f2)) drho
