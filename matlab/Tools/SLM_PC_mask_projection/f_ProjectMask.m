%% Plot Phase Mask either on the PC or on the SLM
function [wrapMask,wrapMaskFig] = f_ProjectMask(r,mask,gl,phaseValues,mingl,...
maxgl,levShft,normMag,binMask,binv,monitorSize,scrnIdx,tit,str, ...
coordType,abs_ang,MaxMask,plotMask) % [wrapMask,wrapMaskFig]
% Inputs:
%  r: polar coordinate
%  mask: function to be plotted. It is wrapped on [-pi,pi] if abs_ang = 2.
%        Complex structure that has not been truncated.
%        mask = exp(i*UnwrappedMask)
%  gl: number of grey levels (normally 256)
%  phaseValues: discretized phi vector on [-pi,pi].
%  mingl,maxgl: minimum/maximum gray level depth. Ref: 0,255
%  levShft: corresponds to the brightness or constant shift of the gl's
%  normMag: normalize magnitude. yes(1); no(0)
%  binMask: binarizes the mask w.r.t the max and min of the phase (boolean)
%  monitorSize: size of the selected screen for coordType = 2 or of the 
%  grid (sSize) for coordType = 1 
%  binv: binary inversion of the mask: yes(1); no(0). Only applies when 
%        binMask=1. It is usefull to be applied for odd p's on LG beams
%  monitorSize: size of the selected screen 
%  scrnIdx: screen number selector. In [1,N] with N the # of screen
%  tit: plot title
%  str: colorbar string when abs_ang = 0
%  coordType: type of calculation of the spatial coordinates. def: 2 
%    -1: size defined by the user, space support defined by the SLM to use
%    -2: size defined by the resolution of the selected screen    
%  abs_ang: custom(0)[str has to be defined for this case], magnitude
%           (1) or phase (2) plot. Doesn't apply for Zernike and LG +
%           Zernike.
%  MaxMask: defines if the mask should be maximized when coordType = 1
%           -0: custom-size mask that depends on the variable sSize   
%           -1: maximizes the mask for coordType = 1
%           -2: maximized mask but keeping its rectangular fashion
%  plotMask:  no (0); on the screen (1); on the SLM (2); on the screen, but
%             a surface (3)
%
% Outputs:
%  wrapMask: truncation and angle operations on mask.
%  wrapMaskFig: figure handler if needed outside the function
%
% Notes:
%  Image is shown with gl gray levels

%% Wrapping and Circular pupil Application
switch abs_ang
  case 0 % No operation, custom input (assumed to be non complex)
    wrapMask = mask;
    figtit = 'Mask';
    if isreal(mask) == false
      warning(['When you choose abs_ang = 0, mask must be selected to be' ...
        'real-valued. Selecting the real part...']);
      wrapMask = real(mask);
    end
    % str: defined in the input
  case 1 % Amplitude
    wrapMask = abs(mask); % Actually, this is an amplitude filter
    figtit = 'Amplitude Mask';
    str = 'Value of amplitude'; % Colorbar string
    %% Normalization constants (amplitude)
    if normMag == 1 % Phase is not changed
      norm = max(wrapMask(:)); % Max value
      wrapMask = wrapMask/norm; % Normalization of the magnitude
    end
  
  case 2 % Phase
   % Circular pupil and wrapping   
   [wrapMask,customMap] = f_MaskWrapCircDiscret(r,mask,binMask,binv,phaseValues,mingl, ...
                                       maxgl,levShft,coordType,plotMask);
   figtit = 'Phase Mask';
   str = 'Wrapped phase value'; % Colorbar string
 end

%% Plot
switch plotMask
 case 0 
  % Won't plot at all
  wrapMaskFig = figure('Visible','off');
 case 1 % PC Screen: normal plot
  % OLD:   
  % slmhfig = figure('color','white','units','normalized','position',...
  % [0 0 1 1],'outerposition',[1/2 0 1/2 1],'Name',tit);
  %% Plot the mask
  wrapMaskFig = figure('color','white','units','normalized','position',...
               [0 0 1 1],'outerposition',[5/10 1/10 1/3 1/2],'Name',figtit); 
  imagesc(wrapMask); axis square; colormap(customMap);
  title(tit);
  set(gca,'xtick',[]); set(gca,'ytick',[]) % No axes values
  cbh = colorbar; cbh.Label.String = str;
  pax = gca; pax.FontSize = 16; % Font size   
  
  %% Add the -pi -> +pi ticks
  if abs_ang == 2
    % Next, pi ticks are added to the colorbar if the phase value extreme
    % points are near -pi and pi; otherwise the default colorbar is shown
    upTol = abs(max(wrapMask(:)) - pi); % upper tolerance
    lowTol = abs(abs(min(wrapMask(:))) - pi); % lower tolerance
    % Tol and a were found empirically
    Tol = 0.3; % Tolerance for pi ticks
    a = 0.15; % Colorbar custom tick adjustment
    if upTol < Tol && lowTol < Tol % Tolerances near pi values
      % cbh.Ticks = linspace(1, 7, 2*pi); % 2*pi: period of e^(i*x) 
      cbh.Ticks = linspace(-pi+a, pi-a, 5);
      cbh.TickLabels = {'-\pi' '-\pi/2' '0'  '\pi/2' '\pi'};
    end 
  end
    
 case 2  % Plot on the SLM. 
  %% Monitor selection and resolution retrieval
  enablechange = true; % SLM figure display monitor activated
  [res,~]= f_changeProjectionMonitor(scrnIdx,enablechange);
  % 2018a: gcf.WindowState = 'maximized';
  % https://blogs.mathworks.com/pick/2018/07/13/maximize-your-figures/
  
  %% Figure handler definitions
  wrapMaskFig = figure('Visible','off','MenuBar','none','Toolbar', ...
                      'none','NumberTitle','off');
  % Hide Menu bar and Tool bar
  % wrapMaskFig.Units = 'Pixels'; % 'color','black', NOT NEEDED FOR NOW
  
  %% Figure size adjusting with a monitor/mask scaling for coordType = 1 or
  %%% a maximized figure for coordType = 2
  if coordType == 1 
    maskSize = monitorSize; % As defined in f_DefineSpace (square-size)
    % res = the real "monitorSize", taken on the "%% Monitor selection and 
    % resolution retrieval" Section   
    switch MaxMask 
       case 0 % Custom-size mask that depends on the variable sSize
         xMov = maskSize(1); % x movement
         yMov = maskSize(2); % y movement
         MidVectMonitor = ceil((res+1)/2); % SLM monitor mid vector
         MidVectMask = ceil((size(wrapMask)+1)/2); % Mask mid vector        
         offsetPixel = MidVectMonitor - MidVectMask; % Pixel position of 
                                                     % the mask
       case 1 % Maximizes the mask in a rectangular screen
         xMov = res(1); % x movement
         yMov = res(2); % y movement
         offsetPixel = [0 0]; % left bottom part
         
       case 2 % Maximizes the mask but keeping its rectangular fashion
         xMov = min(res); % Smallest width of the screen's resolution
         yMov = xMov; % y movement
         MidVectMonitor = [ceil((max(res)+1)/2) 0];
         MidVectMask = [ceil((max(xMov)+1)/2) 0];
         offsetPixel = MidVectMonitor - MidVectMask; % Pixel position of 
                                                   % the mask
    end 
  else % coordType = 2 
    xMov = monitorSize(1); % x movement
    yMov = monitorSize(2); % y movement
    offsetPixel = [1,1]; % Mandatory: pixels have this origin. 
                        % [0,0] doesn't exist
    % full-screen size figures are always created here
  end
  %%% For cordType = 1 and 2:
  set(gca,'Units','Pixels'); % Axis units. get current axis command
  set(gca,'Position',[offsetPixel xMov yMov]); % Figure position
  % [left bottom width height]
  
  %% Figure plotting
  imagesc(wrapMask);  % Plots in SLM screen 
  axis off; colormap(customMap);
  % axis fill;
  wrapMaskFig.Visible = 'on';
  [~] = f_changeProjectionMonitor('Restore'); % Restore default figure
  % OLD
  % set(gca,'xtick',[]); set(gca,'ytick',[]) % No axis values
    
 case 3 % Screen: surface plot
  wrapMaskFig = figure('color','white','units','normalized','position',...
               [0 0 1 1],'outerposition',[5/10 1/10 1/3 1/2],'Name',figtit);
  surf(wrapMask), colormap(gray(gl)), shading interp; % 3D Surface
  axis square; title(tit);
  cbh = colorbar; cbh.Label.String = str;
  pax = gca; pax.FontSize = 16; % Font size
  % No axis values:
  set(gca,'xtick',[]); set(gca,'ytick',[]);  set(gca,'ztick',[]) 
  axis off
end % of switch MaxMask 
end