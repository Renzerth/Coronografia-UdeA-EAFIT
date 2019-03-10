function [energy,radialIntensity] = f_calculateEEF(distribution,shiftCart,metricProfile)
%
% Inputs:
%  distribution: either a simulated 2D image or a camera-taken 2D image
%  shiftCart: [yshift,xshift], works when shiftMask = 1. Percentages of
%             movement of the total size of the mask (cartesian coordinates
%             convention). Ranges per shift: [0,100] (percentage)  
%  metricProfile: 1: vertical profile; 2: horizontal profile
% Outputs:
%

%% Mid and max points of the mask (with the shift application)
% The signs of the shifts account for the cartesian coordinates convention
% [maxX, maxY] = size(distribution);
% maxX = maxX - shiftCart(1); % Shifted X for the SLM 
% maxY = maxY + shiftCart(2); % Shifted Y for the SLM

% midX = round((maxX+1)/2) + mod(maxX,2);
% midY = round((maxY+1)/2) + mod(maxY,2);

%% Profile
[Hprof,Vprof] = f_makeImageProfile(x,y,distribution,tol,shiftCart,tit, ...
                                   plotData,plotH,plotV,oneSideProfile);
                     
%% Profile choosing
switch metricProfile
    case 1 % Vertical profile
        radialIntensity = Vprof;
    case 2 % Horizontal profile
        radialIntensity = Hprof;
    otherwise
        error('metricProfile must be either 1 or 2');
end


% radialIntensity = improfile(distribution,[1,maxY],[midX,midX]);
% radialIntensityOneside = radialIntensity(fix(end/2):end);

%% Enclosed Energy Factor
energy = cumsum(radialIntensityOneside ); % Discrete integration
energy  = energy/energy(end); % Same as normalizing with the max
normIntensity = radialIntensityOneside./max(radialIntensityOneside);

figure('color','white');
plot(1:midX,energy'); hold on
plot(1:midX,normIntensity); hold off
title('Encircled Energy Distribution of Intensity'); grid on;
xlabel('\lambda/D'); ylabel('Relative throughput')
legend({'Encircled Energy Factor', 'Intensity'});
end