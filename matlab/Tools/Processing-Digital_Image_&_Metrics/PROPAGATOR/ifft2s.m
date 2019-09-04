function [signalSynthesis] = ifft2s(signalArray,NM,dxdy)
signalSynthesis = ifftshift(ifft2(fftshift(signalArray))).*dxdy.*NM;
end