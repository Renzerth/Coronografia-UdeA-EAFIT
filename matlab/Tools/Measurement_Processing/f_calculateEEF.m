function [energy] = f_calculateEEF(radialIntensity, cartcoord, titprof, ...
                                               tit,xlab)
% Plots the Enclosed Energy Factor (EEF) along with its respective
% intensity distribution
%
% Inputs:
%  -distribution: either a simulated 2D image or a camera-taken 2D image
%  -shiftCart: [yshift,xshift], works when shiftMask = 1. Percentages of
%             movement of the total size of the mask (cartesian coordinates
%             convention). Ranges per shift: [0,100] (percentage)  
%  metricProfile: 1: vertical profile; 2: horizontal profile
% Outputs:
%  -energy
%  -radialIntensity
%

%% Enclosed Energy Factor (EEF)
energy = cumsum(radialIntensity); % Discrete integration
energy  = energy/energy(end); % Same as normalizing with the max
normIntensity = radialIntensity./max(radialIntensity);

%% Plot of the EEF and its corresponding intensity pattern
figure('color','white');
[hAxes,hLine1,hLine2] = plotyy(cartcoord,energy,cartcoord,normIntensity);
set(hLine1, 'Color','b','LineStyle','-'); set(hLine2, 'Color','r','LineStyle','-');
set(hAxes(1),'YColor','b'); set(hAxes(2),'YColor','r'); set(hAxes,'FontSize',10,'FontWeight','bold');
axis(hAxes(1), 'square'); axis(hAxes(2), 'square'); box(hAxes(1), 'on');
title(strcat(tit,{' '},titprof)); grid on;
xlabel(xlab);
ylabel(hAxes(1),'Relative throughput','FontSize',14,'FontWeight','bold');
ylabel(hAxes(2),'Intensity pattern','FontSize',14,'FontWeight','bold'); grid on;
lgdHandler = legend({'Encircled Energy Factor', 'Intensity'}); lgdHandler.FontSize = 10;

end