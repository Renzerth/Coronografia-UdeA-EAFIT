function [EEFrejection, EEFattenuation] = f_calculateEEFAttenuat(nonCoronagraphicEEF,coronagraphicEEF)
%% Compute 'leakage' of with-vortex and vortex intensities
EEFrejection = nonCoronagraphicEEF./coronagraphicEEF;
EEFattenuation = 1./EEFrejection;
end