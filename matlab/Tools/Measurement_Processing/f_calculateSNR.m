function [SNR] = f_calculateSNR(Measurement, RefMeasurement, cartcoord,...
                                               tit1,tit2,xlab)
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

%% Plot of each intensity in log scale
figure('color','white');
semilogy(cartcoord,Measurement); hold on
semilogy(cartcoord,RefMeasurement); hold off; 
title(tit2); grid on;
xlabel('Angular position [\lambda/D]'); ylabel('Relative contrast of the radial intensities [logscale]')
legend({'Measurement', 'Reference'});

%% SNR
SNR = log10(Measurement) - 0.5*log10(RefMeasurement);

%% Plot of the SNR
figure('color','white');
plot(cartcoord,SNR);
xlabel(xlab); ylabel('SNR')
title(tit1); grid on;

end