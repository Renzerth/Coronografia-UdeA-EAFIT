function [energy,radialIntensity] = f_calculateEEF(distribution,n,PP,M,f,shiftCart,metricProfile,tit)
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

tol = 0; % 0: no need to symmetrically truncate the profile
plotData = 1; % Not needed for the metric
plotH = 0; % Not needed for the metric
plotV = 0; % Not needed for the metric
oneSideProfile = 1; % Specifically needed for this metric

%% Cartesian coordinates
[ySize, xSize] = size(distribution);
xpix = 1:xSize;
ypix = 1:ySize;
xang = f_scalePix2DiffAng(xpix,n,PP,M,f);
yang = f_scalePix2DiffAng(ypix,n,PP,M,f);

% shiftCart = [100,100];
%% Profile
[Hprof,Vprof,~,~,midX,midY] = f_makeImageProfile(xang,yang,distribution,tol,shiftCart,tit, ...
                                   plotData,plotH,plotV,oneSideProfile);
                     
%% Profile choosing
switch metricProfile
    case 1 % Vertical profile
        radialIntensity = Vprof; % One-sided
        midCart = midY;
    case 2 % Horizontal profile
        radialIntensity = Hprof; % One-sided
        midCart = midX;
    otherwise
        error('metricProfile must be either 1 or 2');
end


% radialIntensity = improfile(distribution,[1,maxY],[midX,midX]);
% radialIntensityOneside = radialIntensity(fix(end/2):end);

%% Enclosed Energy Factor
energy = cumsum(radialIntensity); % Discrete integration
energy  = energy/energy(end); % Same as normalizing with the max
normIntensity = radialIntensity./max(radialIntensity);

figure('color','white');
plot(1:midCart,energy'); hold on
plot(1:midCart,normIntensity); hold off
title('Encircled Energy Distribution of Intensity'); grid on;
xlabel('\lambda/D'); ylabel('Relative throughput')
legend({'Encircled Energy Factor', 'Intensity'});
end