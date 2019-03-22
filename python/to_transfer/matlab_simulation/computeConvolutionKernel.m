function [convKernel] = computeConvolutionKernel(freqMeshX,freqMeshY,freqCirc, k, z, planeSize,kernelType)
%% Propagation Kernel - Spectral Transfer functions
switch kernelType
    case 1
        lambda = 2*pi/k;
        convKernel =freqCirc.*exp(1i*k*z*sqrt(1-(lambda*freqMeshX).^2-(lambda*freqMeshY).^2)); % Angular Spectrum -> Direct spectrum Kernel
    case 2
        lambda = 2*pi/k;
        convKernel = exp(1i*k*z)*exp(-1i*pi*lambda*z*(freqMeshX.^2 + freqMeshY.^2)); % Transfer Function -> Direct spectrum Kernel
    case 3
        waveXComponents = 2*pi*(freqMeshX/planeSize(1)); % kx trasverse wave component
        waveYComponents = 2*pi*(freqMeshY/planeSize(2)); % ky trasverse wave component
        waveZSquaredCo = k^2 - waveXComponents.^2 - waveYComponents.^2;
        waveZSquaredCo(waveZSquaredCo < 0) = 0; % Evanescent waves propagration supression
        waveZComponents = sqrt(waveZSquaredCo); % kz axial wave component
        convKernel = exp(1i*waveZComponents*z);
    otherwise
        error('Option is not available.')
end