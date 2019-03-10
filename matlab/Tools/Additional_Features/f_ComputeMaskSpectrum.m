function f_ComputeMaskSpectrum(x,y,mask,maskFTlog,FTmask)
% Computes the spectrum of a compelx mask of the form 
% mask = exp(i*UnwrappedMask), by using the cartesian coordinate vectors
% x,y and by optionally plotting with log scale. FTmask and maskFTlog are
% booleans

if FTmask == 1
%% Inputs
if maskFTlog == 1
 tol = 1; % Tolerance of each side when plotting the FFT profiles since 
          % some values may tend to -infinity when masFTlog = 1
else
 tol = 0; % tol is an index, so it must be an integer
end

 %% Fourier transform
 FFT2D = @(s) ifftshift((fft2(fftshift(s)))); % 2D Fourier Transform
 maskFFT = FFT2D(mask); % FT of the mask (not wrapped)
 switch abs_ang
     case 0
     case 1
     case 2
     maskFFT = abs(maskFFT); % Magnitude of the FT
     if maskFTlog == 1
       logmask = 20*log10(maskFFT); % Magnitude squared of the FT in log scale
       mask = logmask; % Logarithm of the magnitude squared 
       tit = 'Log of the FT of the Mask';
     else
       mask = maskFFT; % Magnitude
       tit = 'FT of the Mask';
     end
 end    
 figure; imagesc(mask); title(tit); % colormap(hot);

 %% Mid and max points of the mask
 [maxX, maxY] = size(mask);
 midX = round((maxX+1)/2) + mod(maxX,2);
 midY = round((maxY+1)/2) + mod(maxY,2);

 %% Profiles in the 2D image
 hold on
 % Horizontal:
 line([1,maxY],[midX,midX],'LineWidth',3,'Color','blue','LineStyle','--'); 
 % Vertical:
 line([midY,midY],[1,maxX],'LineWidth',3,'Color','red','LineStyle','--'); 
 hold off

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
 
  %% OLD
 % a = abs(FFT2D((exp(1*1i*phi))));
 % b = 20*log10(a);
 % figure;
 % imagesc(b);
 % figure;
 % g = improfile(b,[1,1023],[512,512]);
 % plot(g)
end