%% Inputs
tol = 1;

%% Fourier transform
FFT2D = @(s) ifftshift((fft2(fftshift(s)))); % 2D Fourier Transform
maskSpectrum = FFT2D(mask); % FT of the mask (not wrapped)
maskSpectrum = abs(maskSpectrum); % Magnitude of the FT
logmask = 20*log10(maskSpectrum); % Magnitude squared of the FT in log scale
abs_ang_FT = 1; % Magnitude is always plotted
plotType = 1;
f_fig_maskSLM(x,y,r,logmask,gl,glphi,mingl,maxgl,levShft,abs_ang,binMask, ...
              monitorSize,scrnIdx,plotType);

            %% Mid and max points of the mask
[maxX, maxY] = size(logmask);
midX = ceil((maxX+1)/2);
midY = ceil((maxY+1)/2);

%% Profiles in the 2D image
hold on
line([1,maxY],[midX,midX]); % Horizontal
line([midY,midY],[1,maxX]); % Vertical
hold off

%% Horizontal profile
figure;
f = improfile(logmask,[1,maxY],[midX,midX]);
plot(x(1+tol:end-tol),f(1+tol:end-tol)); 
title('Horizontal profile of abs squared FFT');

%% Vertical profile
figure;
g = improfile(logmask,[midY,midY],[1,maxX]);
plot(y(1+tol:end-tol),g(1+tol:end-tol)); 
title('Vertical profile of abs squared FFT');