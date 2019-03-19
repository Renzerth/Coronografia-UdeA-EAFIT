function [energy,Measurement] = f_calculateSNR(x,y,distribution, refdistribution, ...
                                               shiftCart,metricProfile,tit,xlab,ylab,plotData,plotH,plotV,oneSideProfile,dcShift,tol)
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

%% Profile
[x,y,Hprof,Vprof,~,~,~,~] = f_makeImageProfile(x,y,distribution,...
tol,shiftCart,tit,xlab,ylab,plotData,plotH,plotV,oneSideProfile,dcShift);
                     
%% Profile choosing
switch metricProfile
    case 1 % Vertical profile
        Measurement = Vprof; % One-sided
        cartcoord = y;
        titprof = '(vertical profile)';
    case 2 % Horizontal profile
        Measurement = Hprof; % One-sided
        cartcoord =x;
        titprof = '(horizontal profile)';
    otherwise
        error('"metricProfile" must be either 1 or 2');
end

%% SNR

%% Plot of each intensity in log scale
figure('color','white');
semilogy(cartcoord,Measurement); hold on
semilogy(cartcoord,refdistribution); hold off; 
title(strcat(tit,{' '},titprof)); grid on;
xlabel('Angular position [\lambda/D]'); ylabel('Relative contrast (logscale) [Radial intensity]')
legend({'Measurement', 'Reference'});

%% Plot of the SNR

SNR = log(Measurement) - 0.5*log(refdistribution);
plot(cartcoord,SNR);
xlabel('Angular position [\lambda/D]'); ylabel('SNR')
title('SNR')); grid on;
end