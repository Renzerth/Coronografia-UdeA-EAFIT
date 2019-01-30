import pyfftw
import numpy
import time
import scipy
from multiprocessing import cpu_count

#pyfftw.forget_wisdom()
pyfftw.interfaces.cache.disable()

mainDataSize = 536
fastLenght = pyfftw.next_fast_len(mainDataSize)

pyfftw.config.NUM_THREADS = cpu_count()
pyfftw.config.PLANNER_EFFORT = 'FFTW_ESTIMATE'

print(fastLenght)

dataSize = (128,fastLenght,fastLenght)


#f = pyfftw.n_byte_align_empty((127,512,512),16, dtype='complex128')
f = pyfftw.empty_aligned(dataSize, dtype='complex128', n=16)
f[:] = numpy.random.randn(*f.shape)

# first call requires more time for plan creation
# by default, pyfftw use FFTW_MEASURE for the plan creation, which means that many 3D dft are computed so as to choose the fastest algorithm.
fftf=pyfftw.interfaces.numpy_fft.fftn(f)

#help(pyfftw.interfaces)
tas = time.time()
fftf=pyfftw.interfaces.numpy_fft.fftn(f) # here the plan is applied, nothing else.
tas = time.time()-tas
print("3D FFT, pyfftw:", tas)

#f = pyfftw.n_byte_align_empty((127,512,512),16, dtype='complex128')
f = pyfftw.empty_aligned(dataSize, dtype='complex128', n=16)
f[:] = numpy.random.randn(*f.shape)


tas = time.time()
fftf=numpy.fft.fftn(f)
tas = time.time()-tas
print("3D FFT, numpy:", tas)

tas = time.time()
fftf=scipy.fftpack.fftn(f)
tas = time.time()-tas
print("3D FFT, scipy/fftpack:", tas)

# first call requires more time for plan creation
# by default, pyfftw use FFTW_MEASURE for the plan creation, which means that many 3D dft are computed so as to choose the fastest algorithm.
#f = pyfftw.n_byte_align_empty((128,512,512),16, dtype='complex128')
f = pyfftw.empty_aligned(dataSize, dtype='complex128', n=16)
fftf=pyfftw.interfaces.numpy_fft.fftn(f)

tas = time.time()
fftf=pyfftw.interfaces.numpy_fft.fftn(f) # here the plan is applied, nothing else.
tas = time.time()-tas
print("3D padded FFT, pyfftw:", tas)