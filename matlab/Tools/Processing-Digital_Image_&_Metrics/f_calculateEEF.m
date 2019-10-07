function [energy, normIntensity] = f_calculateEEF(radialIntensity, varargin)
% Plots the Enclosed Energy Factor (EEF) along with its respective
% intensity distribution
%
% Inputs:
%  -distribution: either a simulated 2D image or a camera-taken 2D image
%
% Optional Inputs (Enables plotting)
%  -cartcoord: Reference data vector to
%  -titprof: Detail of Title of plot (Profile Direction)
%  -tit: Title of plot
%  -xlab: x axis label

% Outputs:
%  -energy: Radial Enclosed Energy Factor
%  -normIntensity: Normalized Intensity
%% Program Settings and input check
if nargin == 4 && prod(cellfun(@(inputs) ischar(inputs), varargin(2:end))) && isnumeric(varargin{1})
    cartcoord = varargin{1};
    titprof = varargin{2};
    xlab = varargin{3};
    fontSize = 12;
    lineWidth = 1.0;
    plottingEnabled = true;
else
    plottingEnabled = false;
end

%% Enclosed Energy Factor (EEF)
if isnan(radialIntensity(end)); radialIntensity(end) = 0; end; % Remove nans after averaged 0/0 values located at the end
energy = cumsum(radialIntensity); % Discrete integration
energy  = energy/energy(end); % Same as normalizing with the max
normIntensity = radialIntensity./max(radialIntensity);

%% Plot of the EEF and its corresponding intensity pattern
if plottingEnabled
    f_plotEEF(cartcoord, energy, normIntensity, titprof, xlab, fontSize, lineWidth)
end
end