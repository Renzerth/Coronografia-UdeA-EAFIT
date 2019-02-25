function [customMap,discretMap] = f_discretizeMask(phaseValues,phaseMap)
% Inputs:
%  phaseValues:a vector of custom size with phase values ranging on [0,2*Pi]
%  phaseMap: real-valued structure of an phase mask
%
% Outputs:
%  customMap:
%  discretMap:

%% Computes gray levels associated to phase values
grayValues = phaseValues/(2*pi); % Change phase values to a normalized gray map channel
customMap = repmat(grayValues',[1,3]); % Replicates a gray map into RGB channels
NG = length(grayValues); % Computes the quantization level to be used in the phase map

%% Applies modulo operator to phase map within specified phase values
wrappedMap = f_wrapToRange(phaseMap, min(phaseValues), max(phaseValues)); % wraps phase map into the defined span range of phase
discretMap = 2*pi*floor(wrappedMap/(2*pi/NG))/(NG-1) - pi; % Discretization of the phase map trough truncation of 2Pi/NG
% It's span range is on [-Pi,+Pi]
end
