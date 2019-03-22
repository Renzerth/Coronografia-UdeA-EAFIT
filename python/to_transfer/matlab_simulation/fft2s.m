function [signalAnalysis] = fft2s(signalArray,dxdy)
signalAnalysis = fft2(fftshift(signalArray)).*dxdy;
end