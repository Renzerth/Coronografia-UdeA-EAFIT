%% Plot Phase Mask 
function[] = f_fig_maskSLM(x,y,r,mask,m,n,a,b,gl,abs_ang,binMask,plotMask)
% Inputs:
%  x,y: cartesian coordinates vector
%  r: polar coordinate (in cm)
%  mask: function to be plotted. It is wrapped on [-pi,pi] if abs-ang = 2
%  m: y-pos; positive up
%  n: x-pos; positive right
%  a: x-scale
%  b: y-scale
%  gl: number of grey levels (normally 256)
%  abs_ang: Magnitude (1); Phase (2)
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  plotMask:  no (0); on the screen (1); on the SLM (2); on the screen, but
%             a surface (3)
%             plotMask = show; for 0 and 1.
%
% Notes:
%  Image is shown with gl gray levels
%  m,n,a,b only work for plot = 1

%% Wrapping and Circular pupil Application
if abs_ang == 2
    wrappedMask = f_mask_circ_angle_gl(r,mask,binMask);
    tit = 'Phase Mask';
    str = 'Value of phase';
else % abs_ang == 1
    wrappedMask = abs(mask); % Actually, this is an amplitude filter
    tit = 'Amplitude Mask';
    str = 'Value of amplitude';
end

%% Plot
switch plotMask
  case 1 % Screen: normal plot
    figure('color','white','units','normalized','position',...
           [0 0 1 1],'outerposition',[1/2 0 1/2 1],'Name',tit);
    imagesc(x,y,wrappedMask); axis square; colormap(gray(gl));
    title(tit);
    set(gca,'xtick',[]); set(gca,'ytick',[]) % No axes values
    cbh = colorbar; cbh.Label.String = str;
    pax = gca; pax.FontSize = 16; % Font size   
    
  case 2  % Plot on SLM
     hFigure = figure('color','black','units','normalized','position',...
                      [0 0 1 1],'outerposition',...
                      [1.2340+m/100 0.169-n/100 a b]);
                      %[1.2340+m/100 0.169-n/100 0.6*a 0.8*b]);
                         
                    % [xN yN widthN heightN]; N: Normalized
    set(gca,'xtick',[]); set(gca,'ytick',[]) % No axis values
    set(hFigure, 'MenuBar', 'none'); % Hide Menu bar
    set(hFigure, 'ToolBar', 'none'); % Hide Tool bar
    imagesc(x,y,wrappedMask); axis fill; colormap(gray(gl)); % Plots in SLM screen
    set(gca,'xtick',[]); set(gca,'ytick',[]) % no axis values
    
  case 3 % Screen: surface plot
      figure('color','white','units','normalized','position',...
           [0 0 1 1],'outerposition',[1/2 0 1/2 1],'Name',tit);
      surf(x,y,wrappedMask), colormap(gray(gl)), shading interp; % 3D Surface
      axis square; title(tit);
      cbh = colorbar; cbh.Label.String = str;
      set(gca,'xtick',[]); set(gca,'ytick',[]);  set(gca,'ztick',[]) % No axis values
      axis off
  case 0 % Won't plot at all
    
end
end