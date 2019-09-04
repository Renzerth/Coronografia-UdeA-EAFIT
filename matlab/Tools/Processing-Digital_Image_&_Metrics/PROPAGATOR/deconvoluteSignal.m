function [deconvolutedSignal] = deconvoluteSignal(transferFunction, signalDisitribution, normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, halfSize)
%% Functions Padding
paddedSignalSpectrum = padarray(fft2s(signalDisitribution,analysisScaling),halfSize,0,'both');
paddedtransferFunction = padarray(transferFunction,halfSize,0,'both');

%% Convolution
spectralProduct = paddedSignalSpectrum./paddedtransferFunction;
spectralProduct(~isfinite(spectralProduct)) = 0;
deconvolutedSignal = ifft2s(spectralProduct(dataCoordX,dataCoordY), normNMFactor, synthesisScaling);
end