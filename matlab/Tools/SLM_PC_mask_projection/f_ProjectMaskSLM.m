%% Plot Phase Mask 
function slmhfig = f_ProjectMaskSLM(r,mask,gl,glphi,mingl,maxgl, ...
                                    levShft,abs_ang,binMask,monitorSize,...
                                    scrnIdx,plotMask)
% Inputs:
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
  str = 'Wrapped phase value';  
end

%% Plot
switch plotMask
  case 0 
      % Won't plot at all
  case 1 % Screen: normal plot
    % slmhfig = figure('color','white','units','normalized','position',...
           % [0 0 1 1],'outerposition',[1/2 0 1/2 1],'Name',tit);
    slmhfig = figure('color','white','Name',tit); 
    imagesc(wrappedMask); axis square; colormap(gray(gl));
    title(tit);
    set(gca,'xtick',[]); set(gca,'ytick',[]) % No axes values
    cbh = colorbar; cbh.Label.String = str;
    pax = gca; pax.FontSize = 16; % Font size   
    if abs_ang == 2
      % Next, pi ticks are added to the colorbar if the phase value extreme
      % points are near -pi and pi; otherwise the default colorbar is shown
      upTol = abs(max(wrappedMask(:)) - pi); % upper tolerance
      lowTol = abs(abs(min(wrappedMask(:))) - pi); % lower tolerance
      % Tol and a were found empirically
      Tol = 0.3; % Tolerance for pi ticks
      a = 0.15; % Colorbar custom tick adjustment
      if upTol < Tol && lowTol < Tol % Tolerances near pi values
           % cbh.Ticks = linspace(1, 7, 2*pi); % 2*pi: period of e^(i*x) 
           cbh.Ticks = linspace(-pi+a, pi-a, 5);
           cbh.TickLabels = {'-\pi' '-\pi/2' '0'  '\pi/2' '\pi'};
      end 
    end
    
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
    surf(wrappedMask), colormap(gray(gl)), shading interp; % 3D Surface
    axis square; title(tit);
    cbh = colorbar; cbh.Label.String = str;
    % No axis values:
    set(gca,'xtick',[]); set(gca,'ytick',[]);  set(gca,'ztick',[]) 
    axis off
  
end
end