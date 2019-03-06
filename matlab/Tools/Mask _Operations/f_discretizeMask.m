function [discretMap,customMap] = f_discretizeMask(phaseValues,phaseMap)
% Here, the discretization range and the phase map's range are equalled
% The phase is mapped in a bounded gray-scale
%
% Inputs:
%  phaseValues:a vector of custom size with phase values ranging on [0,2*Pi]
%  phaseMap: real-valued structure of a phase mask
%
% Outputs:
%  discretMap: discretized version of phase Map ranging on [-Pi,Pi]
%  customMap: map of the generated gray levels in the RGB channels
%
% Created by Juan Jose Cadavid on February 2019


%% Computes the gray levels associated with the phase values
grayValues = phaseValues/(2*pi); % Changes the phase values to a normalized
                                 % gray map channel
customMap = repmat(grayValues',[1,3]); % Replicates a gray map into RGB 
                                       % channels
GL = length(grayValues); % Computes the quantization level to be used in 
                         % the phase map

%% Applies the modulo operator to phase map within  specified phase values
% wraps phase map into the defined span range of phase
wrappedMap = f_wrapToRange(phaseMap, min(phaseValues), max(phaseValues)); 
% Discretization of the phase map trough truncation of 2Pi/NG
discretMap = 2*pi*floor(wrappedMap/(2*pi/GL))/(GL-1) - pi; 
% It's span range is on [-Pi,+Pi] thanks to the "-pi term"
% floor(wrappedMap/(2*pi/GL)): discretization
% (GL-1): normalization factor
end
