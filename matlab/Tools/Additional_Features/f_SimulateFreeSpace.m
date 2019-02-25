function f_SimulateFreeSpace(x,y,X,Y,r,mask,starAmplitude,planetAmplitude, ...
pixelSeparation,w1,w2,rPupilSize,showIin,showPupilin,showFPmag,logscale,...
showFPphas,showPhasout,showMagout,showIout)
%% Vortex Simulation With Gaussians - Coronagraphy
% Inputs:
% x,y: cartesian coordinate vectors
% X,Y,r: meshgrids of cartesian and radial polar
% mask: complex structure that has not been truncated and is 
%  wrapped on [-pi,pi]. mask = exp(i*UnwrappedMask).
% star-planet_Amplitude: of each body. The star is centered and the body
% can be separated by a:
% pixelSeparation: of the two bodies
% w1,w2: with of each body (star-planet)
% rPupilSize:  Input pupil radius. percentage [0,1]; ref: 0.5
% showIin: show input intensity
% showPupilin: show input pupil
% showFPmag: show Fourier plane magnitude
% logscale: apply logscale to the Fourier plane magnitude
% showFPphas: show Fourier plane phase
% showPhasout: show output phase
% showMagout: show output magnitude
% showIout: sow output intensity
%
% Outputs:
%
%%% Part 1, only the input field is generated
% Purpose:
% A vortex is modulated with two bodies with variable size, intensity and
% relative separation.
% The star is brighter, bigger and the vortex is centered on it. The planet
% is the one that one wants to recover with coronagraphy. The filtering is 
% applied on the spatial frequency plane
%
% Notes:
% Star: primary body (one); Planet: secondary body (two)

%% Bodies' relative intensity 
contrast = starAmplitude/planetAmplitude; % Intensities ratios
disp(['Relative contrast between bodies: ' num2str(contrast)]);

%% Bodies' position
disp(['Linear distance between bodies: ' num2str(pixelSeparation)]);
beamCenter1 = [0,0]; % Primary star centered
beamCenter2 = [0,pixelSeparation]; % X-deviation of the secondary body

%% Input bodies' gaussian amplitudes
gBeam = exp(-((X - beamCenter1(1)).^2 + (Y - beamCenter1(2)).^2)/w1^2);
gBeam2 = exp(-((X - beamCenter2(1)).^2 + (Y - beamCenter2(2)).^2)/w2^2);

%% Input Field creation (superposition of bodies)
field =  (starAmplitude*gBeam + planetAmplitude*gBeam2); % Superposition  

%% Input field plot
if showIin == 1
    figure, set(gcf,'color','w');
    imagesc(x,y,abs(field).^2), title('Object plane input intensity'); 
    xlabel('x-axis'), ylabel('y-axis'), cbh = colorbar; % colormap winter; 
    cbh.Label.String = 'Intensity'; % Normalized intensity
    pax = gca; pax.FontSize = 14; axis square;% Font size
end

%% Input pupil of the entrance field
circularFunction = double(r <= rPupilSize);

%% Input pupil plot
if showPupilin == 1
 figure; set(gcf,'color','w');
 imagesc(x,y,circularFunction), title('Pupil apperture'); % Pupil apperture
 xlabel('x-axis'), ylabel('y-axis'), cbh = colorbar; % colormap winter; 
 cbh.Label.String = 'Binary value'; colormap(gray(2)); axis square;
 cbh.Ticks = linspace(0, 1, 2); cbh.TickLabels = {'0' '1'}; % Binary 
 pax = gca; pax.FontSize = 14; % Font size
end

%% 2D FFT definition
FFT2D = @(s) ifftshift((fft2(fftshift(s)))); % 2D Fourier Transform
% Explanation:
%       g0 = fftshift(g); % Shift
%       G0 = fft2(g0)*dA; % 2D fft and dxdy scaling (dx = dy)
%       G = fftshift(G0); % Center
iFFT2D = @(s) ifftshift(ifft2(fftshift(s))); % 2D Inverse Fourier Transform

%% Fourier plane modulation
propagationField = circularFunction.*field; % Entrance pupil field
spectralFiltering = FFT2D(propagationField).*mask; % Fourier plane
                                                   % phase modulation
response = iFFT2D(spectralFiltering); % Propagation through a lens
intensities = abs(response).^2; % Field squared is the intensity

%% Plot Fourier plane magnitude

if showFPmag == 1
 if logscale == 1
  figure, set(gcf,'color','w'); imagesc(x,y,log(abs(spectralFiltering)));
 else
  figure, set(gcf,'color','w'); imagesc(x,y,abs(spectralFiltering));
 end
  title('Magnitude on the Fourier plane (logarithmic scale)'); 
  xlabel('x-axis'), ylabel('y-axis'); cbh = colorbar; % colormap winter; 
  cbh.Label.String = 'Magnitude'; axis square;
  pax = gca; pax.FontSize = 14; % Font size  
end

%% Plot Fourier plane phase
if showFPphas == 1
    figure, set(gcf,'color','w'); imagesc(x,y,angle(spectralFiltering));
    title('Phase on the Fourier plane'); 
    xlabel('x-axis'), ylabel('y-axis'); cbh = colorbar; % colormap winter; 
    cbh.Label.String = 'Phase'; axis square;
    pax = gca; pax.FontSize = 14; % Font size  
end

%% Plot output phase
if showPhasout == 1
    figure, set(gcf,'color','w'); imagesc(x,y,angle(response));
    title('Phase response on the image plane'); 
    xlabel('x-axis'), ylabel('y-axis'); cbh = colorbar; % colormap winter; 
    cbh.Label.String = 'Phase'; axis square;
    pax = gca; pax.FontSize = 14; % Font size  
end

%% Plot output magnitude
if showMagout == 1
    figure, set(gcf,'color','w'); imagesc(x,y,abs(response));
    title('Magnitude response on the image plane'); 
    xlabel('x-axis'), ylabel('y-axis'); cbh = colorbar; % colormap winter; 
    cbh.Label.String = 'Magnitude'; axis square;
    pax = gca; pax.FontSize = 14; % Font size  
end

%% Plot output intensity
if showIout == 1
    figure, set(gcf,'color','w'); imagesc(x,y,intensities);
    title('Intensity response on the image plane'); 
    xlabel('x-axis'), ylabel('y-axis'); cbh = colorbar; % colormap winter; 
    cbh.Label.String = 'Intensity'; axis square;
    pax = gca; pax.FontSize = 14; % Font size  
end