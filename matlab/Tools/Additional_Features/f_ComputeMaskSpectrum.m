function f_ComputeMaskSpectrum(x,y,mask,abs_ang,maskFTlog,FTmask)
% Computes the spectrum of a complex mask with severak features as
% profiles, log scale and phase/angle plots
% 
% Inputs:
% x,y: cartesian coordinate vectors
% mask: function to be plotted. It is wrapped on [-pi,pi] if abs_ang = 2.
%       Complex structure that has not been truncated.
%       mask = exp(i*UnwrappedMask)
% abs_ang: custom(0)[str has to be defined for this case], magnitude
%           (1) or phase (2) plot. Doesn't apply for Zernike and LG +
%           Zernike.
% maskFTlog: boolean for optionally plotting with log scale. 
% FTmask: boolean for activating this function


if FTmask == 1
%% Inputs
if maskFTlog == 1 && abs_ang ~= 0
 tol = 1; % Tolerance of each side when plotting the FFT profiles since 
          % some values may tend to -infinity when masFTlog = 1
else
 tol = 0; % tol is an index, so it must be an integer and is not needed for
          % a non-log scale
end

 %% Fourier transform
 FFT2D = @(s) ifftshift((fft2(fftshift(s)))); % 2D Fourier Transform
 maskFFT = FFT2D(mask); % FT of the mask (not wrapped)
 switch abs_ang
     case 0 % No operation, custom input (assumed to be non complex)  
         warning(['abs_ang = 0 is not valid for FTmask = 1, plotting ' ... 
                 'anyways its real part...']);
         mask = real(maskFFT); % Real part of the FT
         tit = 'Real part of the FT of the Mask';
     case 1 % Amplitude
         mask = abs(maskFFT); % Magnitude of the FT
         if maskFTlog == 1
           logmask = 20*log10(mask); % Magnitude squared of the FT in
                                     % log scale
           mask = logmask; % Logarithm of the magnitude squared 
           tit = 'Log of the FT of the Mask';
         else
           tit = 'Magnitude of the FT of the Mask';
         end
     case 2 % Phase
         mask = angle(maskFFT); % Phase of the FT
         tit = 'Phase of the FT of the Mask';
 end
 
%% Plot the mask and take a horizontal and a vertical profile
plotData = 1; % 1: plots the mask
plotH = 1; % 1: plot the horizontal profile
plotV = 1; % 1: plot the vertical profile
oneSideProfile = 0; % 0: two-sided profile
shiftCart = [0,0]; % No shift here
f_makeImageProfile(x,y,mask,tol,shiftCart,tit,plotData,plotH,plotV, ...
                                oneSideProfile);
end