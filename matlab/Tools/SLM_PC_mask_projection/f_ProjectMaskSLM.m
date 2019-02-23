%% Plot Phase Mask 
function slmhfig = f_ProjectMaskSLM(x,y,r,mask,gl,glphi,mingl,maxgl,levShft,abs_ang,binMask,monitorSize,scrnIdx,plotMask)
% Inputs:
%  x,y: cartesian coordinates vector
%  r: polar coordinate (in cm)
%  mask: function to be plotted. It is wrapped on [-pi,pi] if abs-ang = 2
%  gl: number of grey levels (normally 256)
%  glphi: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  abs_ang: Magnitude (1); Phase (2)
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  monitorSize: size of the selected screen 
%  screenIndex: screen number selector. In [1,N] with N the # of screen
%  plotMask:  no (0); on the screen (1); on the SLM (2); on the screen, but
%             a surface (3)
%             plotMask = show; for 0 and 1.
%
% Output:
%  slmhfig: figure handler if needed outside the function
%
% Notes:
%  Image is shown with gl gray levels
%  m,n,a,b only work for plot = 1

%% Wrapping and Circular pupil Application
switch abs_ang 
 case 0 % No operation, custom input (assumed to be non complex)
  wrappedMask = mask;
  tit = 'Mask';
  str = 'Amplitude';
 case 1 % Amplitude
  wrappedMask = abs(mask); % Actually, this is an amplitude filter
  tit = 'Amplitude Mask';
  str = 'Value of amplitude';
  
 case 2 % Phase
  wrappedMask = f_MaskWrapCircDiscret(r,mask,binMask,glphi,mingl,maxgl, ...
                                    levShft);
  tit = 'Phase Mask';
  str = 'Value of phase';  
end

%% Plot
switch plotMask
  case 1 % Screen: normal plot
    slmhfig = figure('color','white','units','normalized','position',...
           [0 0 1 1],'outerposition',[1/2 0 1/2 1],'Name',tit);
    imagesc(x,y,wrappedMask); axis square; colormap(gray(gl));
    title(tit);
    set(gca,'xtick',[]); set(gca,'ytick',[]) % No axes values
    cbh = colorbar; cbh.Label.String = str;
    pax = gca; pax.FontSize = 16; % Font size   
    
  case 2  % Plot on the SLM
    enablechange = true; % SLM figure display monitor activated
    f_changeProjectionMonitor(scrnIdx,enablechange); % Allow full-screen 
                                                     % size figures
    offsetPixel = [1,1]; % Mandatory: pixels have this origin [0,0] doesn't
                         % exist
    slmhfig = figure('Visible','off','MenuBar','none','Toolbar','none', ...
                     'NumberTitle','off');
    % Hide Menu bar and Tool bar
    slmhfig.Units = 'Pixels'; % 'color','black',
    set(gca,'Units','Pixels');
    set(gca,'Position',[offsetPixel monitorSize(1) monitorSize(2)]);
    %set(gca,'xtick',[]); set(gca,'ytick',[]) % No axis values
    image(wrappedMask);  % Plots in SLM screen 
    axis off; colormap(gray(gl));
    % axis fill;
    slmhfig.Visible = 'on';
    [~] = f_changeProjectionMonitor('Restore'); % Restore default figure
    
  case 3 % Screen: surface plot
    slmhfig = figure('color','white','units','normalized','position',...
                     [0 0 1 1],'outerposition',[1/2 0 1/2 1],'Name',tit);
    surf(x,y,wrappedMask), colormap(gray(gl)), shading interp; % 3D Surface
    axis square; title(tit);
    cbh = colorbar; cbh.Label.String = str;
    % No axis values:
    set(gca,'xtick',[]); set(gca,'ytick',[]);  set(gca,'ztick',[]) 
    axis off
  case 0 % Won't plot at all
    
end
end