import numpy as np

samp = 0.18
lamb = 5e-6
ngrid = 1024
dz = 1000000

xrhosqr = np.tile(((np.arange(ngrid, dtype = np.float64) - int(ngrid/2)) / (ngrid * samp))**2, (ngrid, 1))
rhosqr = xrhosqr + np.transpose(xrhosqr)
rhosqr = np.roll(np.roll(rhosqr, int(ngrid/2), 0), int(ngrid/2), 1)

aa = np.fft.fftshift(np.exp((1j*np.pi*lamb*dz)*rhosqr))

plt.imshow(np.angle(aa))