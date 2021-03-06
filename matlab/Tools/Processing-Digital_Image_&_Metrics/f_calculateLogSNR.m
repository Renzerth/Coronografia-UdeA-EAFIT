function [SNR] = f_calculateLogSNR(Measurement, RefMeasurement, varargin)
% Plots the SNR between a signal and its reference
%
% Inputs:
%  distribution: either a simulated 2D image or a camera-taken 2D image
%  shiftCart: [yshift,xshift], works when shiftMask = 1. Percentages of
%             movement of the total size of the mask (cartesian coordinates
%             convention). Ranges per shift: [0,100] (percentage)
%  metricProfile: 1: vertical profile; 2: horizontal profile
% Outputs:
%  energy
%  radialIntensity

%% Program Settings and input check
if nargin == 5 && prod(cellfun(@(inputs) ischar(inputs), varargin(2:end))) && isnumeric(varargin{1})
    cartcoord = varargin{1};
    titprof = varargin{2};
    xlab = varargin{3};
    fontSize = 12;
    lineWidth = 1.0;
    
    plottingEnabled = true;
else
    plottingEnabled = false;
end

%% Logarithmic SNR
Measurement(Measurement == 0) = 1e-3/255; % remove zeros by a 1000th of a single value over the overall saturated range
SNR = log10(Measurement) - 0.5*log10(RefMeasurement);

%% Plotting
if plottingEnabled
    f_plotLogSNR(cartcoord,RefMeasurement,Measurement,SNR,titprof,xlab,fontSize,lineWidth)
end
end