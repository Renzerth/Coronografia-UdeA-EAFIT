%% Plot Phase Msk on the screen, not the SLM
function pcfig = f_fig_maskPCscreen(x, y, wrappedMask, tit, gl, showM)
% Inputs:
%  x,y: cartesian coordinates vector
%  WrappedMask: function to be plotted. The mask should already be wrapped
%  with the angle function
%  tit: Figure title
%  gl: number of grey levels (normally 256)
%  showM:  no (0); on the screen (1); on the SLM (2)
%
% Output:
%  pcfig: figure handler if needed outside the function
%
% Notes:
%  Image is shown with gl gray levels

if showM == 1
  pcfig = figure('color','white','Name','Phase Mask'); 
  % imagesc(x,y,wrappedMask); axis square; colormap(hot(gl)); % Hot
  imagesc(x,y,wrappedMask); axis square; colormap(gray(gl));% Gray
  title(tit);
  cbh = colorbar; cbh.Label.String = 'Wrapped phase value';
  
  % Next, pi ticks are added to the colorbar if the phase value extreme
  % points are near -pi and pi; otherwise the default colorbar is shown
  upTol = abs(max(wrappedMask(:))-pi); % upper tolerance
  lowTol = abs(abs(min(wrappedMask(:)))-pi); % lower tolerance
  Tol = 0.3; % Tolerance
  a = 0.15; % Colorbar custom tick adjustment
  if upTol < Tol && lowTol < Tol
       % cbh.Ticks = linspace(1, 7, 2*pi); % 2*pi: period of e^(i*x) 
       cbh.Ticks = linspace(-pi+a, pi-a, 5);
       cbh.TickLabels = {'-\pi' '-\pi/2' '0'  '\pi/2' '\pi'};
       set(gca,'xtick',[]); set(gca,'ytick',[]) % No axes values
       pax = gca; pax.FontSize = 14; % Font size
  end 
end
end