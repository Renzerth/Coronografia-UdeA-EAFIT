function [signalSynthesis] = ifft2s(signalArray,NM,dxdy)
signalSynthesis = fftshift(ifft2(signalArray)).*dxdy.*NM;
end