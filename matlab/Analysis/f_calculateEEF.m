function [energy,radialIntensity] = f_calculateEEF(x,y,distribution, ...
                                               shiftCart,metricProfile,tit)
% Plots the Enclosed Energy Factor (EEF) along with its respective
% intensity distribution
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
%
%% Metric-specific default parameters
tol = 0; % 0: no need to symmetrically truncate the profile. Ref: 0
plotData = 1; % Shows the profile lines. Ref: 1
plotH = 0; % Not needed for the metric. Ref: 0
plotV = 0; % Not needed for the metric. Ref: 0
oneSideProfile = 2; % Specifically needed for this metric. Ref: 1
dcShift = 0; % Only used for spectra (Fourier analysis)
xlab = 'Angular position [\lambda/D]';
ylab = 'Angular position [\lambda/D]';

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

%% Plot of the
figure('color','white');
plot(cartcoord,energy); hold on
plot(cartcoord,normIntensity); hold off; 
title(strcat(tit,{' '},titprof)); grid on;
xlabel('Angular position [\lambda/D]'); ylabel('Relative throughput')
legend({'Encircled Energy Factor', 'Intensity'});

end