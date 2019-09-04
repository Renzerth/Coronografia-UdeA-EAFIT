function [convolutedSignal] = convoluteSignal(transferFunction, signalDisitribution, normNMFactor, analysisScaling, synthesisScaling, dataCoordX, dataCoordY, halfSize)
%% Functions Padding
paddedSignalSpectrum = padarray(fft2s(signalDisitribution,analysisScaling),halfSize,0,'both');
paddedtransferFunction = padarray(transferFunction,halfSize,0,'both');

%% Convolution
spectralProduct = paddedtransferFunction.*paddedSignalSpectrum;
convolutedSignal = ifft2s(spectralProduct(dataCoordX,dataCoordY), normNMFactor, synthesisScaling);
end