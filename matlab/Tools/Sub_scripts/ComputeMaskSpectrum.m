%% OLD
% a = abs(FFT2D((exp(1*1i*phi))));
% b = 20*log10(a);
% figure;
% imagesc(b);
% figure;
% g = improfile(b,[1,1023],[512,512]);
% plot(g)

%% Inputs
tol = 1;

%% Fourier transform
FFT2D = @(s) ifftshift((fft2(fftshift(s)))); % 2D Fourier Transform
maskFFT = FFT2D(mask); % FT of the mask (not wrapped)
maskFFT = abs(maskFFT); % Magnitude of the FT
if maskFTlog == 1
    logmask = 20*log10(maskFFT); % Magnitude squared of the FT in log scale
    mask = logmask; % Logarithm of the magnitude squared 
    tit = 'Log of the FT of the Mask';
else
    mask = maskFFT; % Magnitude
    tit = 'FT of the Mask';
end
figure; imagesc(mask); title(tit);
% colormap(hot);

%% Mid and max points of the mask
[maxX, maxY] = size(mask);
midX = ceil((maxX+1)/2);
midY = ceil((maxY+1)/2);

%% Profiles in the 2D image
% hold on
% % Horizontal:
% line([1,maxY],[midX,midX],'LineWidth',3,'Color','blue','LineStyle','--'); 
% % Vertical:
% line([midY,midY],[1,maxX],'LineWidth',3,'Color','red','LineStyle','--'); 
% hold off

%% Horizontal profile
figure;
f = improfile(mask,[1,maxY],[midX,midX]);
plot(x(1+tol:end-tol),f(1+tol:end-tol)); 
title('Horizontal profile of abs squared FFT');

%% Vertical profile
figure;
g = improfile(mask,[midY,midY],[1,maxX]);
plot(y(1+tol:end-tol),g(1+tol:end-tol)); 
title('Vertical profile of abs squared FFT');