%% Fourier transform
FFT2D = @(s) ifftshift((fft2(fftshift(s)))); % 2D Fourier Transform
mask = FFT2D(mask); % FT of the mask (not wrapped)
mask = abs(mask); % Magnitude of the FT
mask = 20*log(mask.^2); % Magnitude squared of the FT in log scale
abs_ang_FT = 1; % Magnitude is always plotted
f_fig_maskSLM(x,y,r,mask,m,n,a,b,gl,abs_ang_FT,binMask,plotMask); % Plot FT

%% Mid and max points of the mask
[maxX, maxY] = size(mask);
midX = ceil((maxX+1)/2);
midY = ceil((maxY+1)/2);  

%% Horizontal profile
figure;
f = improfile(mask,[1,maxX],[midY,midX]);
plot(x,f); 
title('Horizontal profile of abs squared FFT');

%% Vertical profile
figure;
g = improfile(mask,[midY,midX],[1,maxX]);
plot(x,g); 
title('Vertical profile of abs squared FFT');