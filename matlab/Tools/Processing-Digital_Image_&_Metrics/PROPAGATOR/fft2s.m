function [signalAnalysis] = fft2s(signalArray,dxdy)
signalAnalysis = ifftshift(fft2(fftshift(signalArray))).*dxdy;
end