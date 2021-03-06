import pyfftw
import numpy
import time
import scipy
from multiprocessing import cpu_count

#pyfftw.forget_wisdom()
pyfftw.interfaces.cache.disable()

mainDataSize = 1024
volumeSize = 20

fastLenghtT = pyfftw.next_fast_len(mainDataSize)
fastLenghtV = pyfftw.next_fast_len(volumeSize)

pyfftw.config.NUM_THREADS = cpu_count()
pyfftw.config.PLANNER_EFFORT = 'FFTW_ESTIMATE'

print(fastLenghtT)
print(fastLenghtV)

dataSize = (fastLenghtV,fastLenghtT,fastLenghtT)

f = pyfftw.empty_aligned(dataSize, dtype='complex128', n=16)
f[:] = numpy.random.randn(*f.shape) + 1j*numpy.random.randn(*f.shape)

tas = time.time()
fftf=pyfftw.interfaces.numpy_fft.fftn(f) # here the plan is applied, nothing else.
tas = time.time()-tas
print("3D FFT, pyfftw:", tas)

f = pyfftw.empty_aligned(dataSize, dtype='complex128', n=16)

tas = time.time()
fftf=numpy.fft.fftn(f)
tas = time.time()-tas
print("3D FFT, numpy:", tas)

tas = time.time()
fftf=scipy.fftpack.fftn(f)
tas = time.time()-tas
print("3D FFT, scipy/fftpack:", tas)

# https://stackoverflow.com/questions/25812255/row-and-column-headers-in-matplotlibs-subplots
#TCIndex = int(index % TCSize)
#grayIndex = int((index - TCIndex)/TCSize % GLSize)
#ax.set_title('TC: %1.1f::GL: %1d' % (TCRanges[TCIndex],GLRanges[grayIndex]),fontsize=14,position=(0.5,1.0))
#https://stackoverflow.com/questions/40777772/python-integer-division-operator-vs-math-floor