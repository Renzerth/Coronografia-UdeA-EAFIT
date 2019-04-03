function [energy,radialIntensity,cartcoord,titprof] = f_calculateEEF(x,y,distribution, ...
                                               shiftCart,metricProfile,tit,xlab,ylab,plotData,plotH,plotV,oneSideProfile,dcShift,tol)
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

%% Profile
[x,y,Hprof,Vprof,~,~,~,~] = f_makeImageProfile(x,y,distribution,...
tol,shiftCart,tit,xlab,ylab,plotData,plotH,plotV,oneSideProfile,dcShift);
                     
%% Profile choosing
switch metricProfile
    case 1 % Vertical profile
        radialIntensity = Vprof; % One-sided
        cartcoord = y;
        titprof = '(vertical profile)';
    case 2 % Horizontal profile
        radialIntensity = Hprof; % One-sided
        cartcoord =x;
        titprof = '(horizontal profile)';
    otherwise
        error('"metricProfile" must be either 1 or 2');
end

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
xlabel('Angular position [\lambda/D]');
ylabel(hAxes(1),'Relative throughput','FontSize',14,'FontWeight','bold');
ylabel(hAxes(2),'Intensity pattern','FontSize',14,'FontWeight','bold'); grid on;
lgdHandler = legend({'Encircled Energy Factor', 'Intensity'}); lgdHandler.FontSize = 10;

end