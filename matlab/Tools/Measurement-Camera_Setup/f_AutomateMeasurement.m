function [refmeasfullpath] = f_AutomateMeasurement(src,vid,Xslm,Yslm,rSLM, ...
phiSLM,Xpc,Ypc,rPC,phiPC,s,ph0,p,WsizeRatio,L,f_FR,bcst,period,T0, ...
frkTyp,Aalpha,Angalp,Angbet,z_coeff,z_a,z_pupil,z_disp_wrap,z_plot, ...
normMag,binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx,coordType, ...
abs_ang,MaxMask,maskSel,ltcvect,lglvect,totalImgs,waitbeforemeas, ...
imgpartPath,dataformat,imgformat,measfullpath,infoDelim,cameraPlane, ...
tcvect,glvect,measSimulated,recordingDelay,whiteblackref,beepSound)
% Plots phase masks on the Fourier plane of the vortex coronagraph and
% takes images of either its Lyot or PSF plane
%
% Inputs:
%  vid: video input object
%
% Ouputs:
%

%% Automated measurements
expImgs = cell(1,totalImgs); % Cell with the experimental images
MeasInfo = cell(1,totalImgs); % Same initialization as expImgs
idxgral = 1; % Initialization of the general index that runs 
             % on the range: [1,totalImgs]

%% Measurements initialization
% Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
pause(waitbeforemeas); % Seconds before measuring as a safety measurement
t1_dt = datetime; % store time
disp('Measurement started'); disp(t1_dt);

if measSimulated == 0
    initialExposure =  src.Exposure; % Store the current camera's exposure
else
    initialExposure = []; % Empty, won't be used at all
end


%% Measurements
for idxtc = 1:ltcvect 
  
  %% Read a value of the tc vector  
  tc = tcvect(idxtc); % Specific tc for this iteration
  
  %% Dynamically change the camera's exposure with a characterization

  for idxgl = 1:lglvect
    %% Read a value of the gl vector and adjust exposure if needed
     phaseValues = glvect(idxgl); % Specific phase values for this iteration
      
   %% Generate the phase mask and display it on the SLM
    phaseValues = linspace(0,2*pi,phaseValues);
    plotMask = 2; % Select SLM for plotting
    [X,Y,r,phi] = f_SelectCoordinates(Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,...
                                      phiPC,plotMask);
    % OLD: [~,wrapMaskSlm,slmfig,~]                              
    [~,~,slmfig,~] = f_PlotSelectedMask(X,Y,r,phi,phaseValues,...
    tc,s,ph0,p,WsizeRatio,L,f_FR,bcst,period,T0,frkTyp,Aalpha,Angalp, ...
    Angbet,z_coeff,z_a,z_pupil,z_disp_wrap,z_plot,normMag,binMask,binv, ...
    MaskPupil,rSize,monitorSize,scrnIdx,coordType,abs_ang,MaxMask, ...
    plotMask,maskSel);      
    
    %% Display the mask on the PC
    plotMask = 1; % Select PC
    [X,Y,r,phi] = f_SelectCoordinates(Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,...
                                      phiPC,plotMask);                             
    [~,wrapMaskPC,pcfig,~] = f_PlotSelectedMask(X,Y,r,phi,phaseValues,tc,s,ph0,p,...
    WsizeRatio,L,f_FR,bcst,period,T0,frkTyp,Aalpha,Angalp,Angbet, ...
    z_coeff,z_a,z_pupil,z_disp_wrap,z_plot,normMag,binMask,binv, ...
    MaskPupil,rSize,monitorSize,scrnIdx,coordType,abs_ang,MaxMask, ...
    plotMask,maskSel);     
    set(pcfig,'units','normalized','position',[6/11 2/10 3/7 1/2]);
    
    %% Wait between cycles
    pause(recordingDelay/2); % Displays the mask for these seconds 
  
   %% Record a snapshot
    if measSimulated == 0
        snap = getsnapshot(vid); % Real measurements
        wait(vid); % Waits until vid is not running or logging
    else % measSimulated = 1
          % OLD: snap = wrapMaskSlm; 
          snap = wrapMaskPC; % the mask is saved
    end
    
   %% Displaying the camera on the PC
    % The numbers after 'position' were empirically obtained
    camfig = figure('units','normalized','position',[1/12 2/10 3/7 1/2]);
                                               % [left bottom width height]
    imagesc(snap); % normalized
    colorbar; title(strcat('Camera image:',{' '},cameraPlane));
    pax = gca; pax.FontSize = 16; % Font size
    % The numbers after 'position' were empirically obtained
    
    %% Time for the snapshot display on the pc
    pause(recordingDelay/2); % seconds for the image to be shown on the pc
    
   %% Saving the measurement
    if measSimulated == 0
         expImgs{idxgral} = im2double(snap); % An extructure with all the images
    else % measSimulated = 1
         % If you want to simulate with the shifted mask, put wrapMaskslm
         expImgs{idxgral} = snap; % Saves the mask: no double conversion 
         % needed since it has the phase values
    end 
    
    tcstr = strcat('tc',infoDelim,num2str(tcvect(idxtc))); 
    glstr = strcat('gl',infoDelim,num2str(glvect(idxgl)));
    MeasInfo{idxgral} = strcat(tcstr,infoDelim,glstr); % Dataname for each 
                                           % experimental data
    imgfullpath = strcat(imgpartPath,MeasInfo{idxgral});
    if measSimulated == 0
        % Explanation: imwrite(variables,directory+filename+extension)
        imwrite(snap, strcat(imgfullpath,dataformat)); % snap is uint8
    else % measSimulated == 1
        % Explanation: saveas(variable,directory+filename,extension)
        saveas(pcfig,imgfullpath,imgformat); % Saves the last shown figure
        % Other options: savefig and print.
    end
     
    %% Preparation for a new measurement iteration          
    stridxgral = num2str(idxgral); strtotalImgs = num2str(totalImgs);
    disp(strcat(stridxgral,' out of ',{' '},strtotalImgs, ...
         ' images recorded'));
    idxgral = idxgral + 1; % The general index increases   
    
    % Close the shown figures  
    if ishandle(pcfig)
      close(pcfig);
    end
    if ishandle(camfig)
      close(camfig); 
    end
    if ishandle(slmfig{1})
      close(slmfig{1});
    end
  end
end

%% Store all measurements in a .mat file
% Explanation: save(directory+filename,variables) % ,'-append'
save(measfullpath,'expImgs','MeasInfo'); % Save as .mat. expImgs is saved as double

%% Reference measurement
% For tc=0 and for a gray level of 0 (black) or 256 (white)
% This is the non-coronagraphic PSF reference: the system response without
% a phase mask

%%% Generate reference mask:
tcref = 0; % Always null: no OAM
if measSimulated == 0 % Measurement
    whiteblackref = 0;
    plotMaskref = 2; % Select SLM for plotting
else % Simulation
    whiteblackref = 1;
    plotMaskref = 1; % Select PC for plotting
end

switch whiteblackref
    case 0
         glref = 0;
         phaseValuesref = [0 1]; % Black
    case 1
         glref = 256;
         phaseValuesref = [0 255]; % White 
end
phaseValuesref = phaseValuesref*2*pi/255; % Conversion from gray-levels to 
                                 % phase levels      
pref = 0; % 0: No radial nodes
binMaskref = 1; % 1: mask is binarized
maskSelref = 1; % LG mask
WsizeRatioref = 100; % Radial nodes' width of 100
[Xref,Yref,rref,phiref] = f_SelectCoordinates(Xslm,Yslm,rSLM,phiSLM, ...
                                            Xpc,Ypc,rPC,phiPC,plotMaskref);
[~,refMask,reffig,~] = f_PlotSelectedMask(Xref,Yref,rref, ... 
phiref,phaseValuesref,tcref,s,ph0,pref,WsizeRatioref,L,f_FR,bcst, ...
period,T0,frkTyp,Aalpha,Angalp,Angbet,z_coeff,z_a,z_pupil, ...
z_disp_wrap,z_plot,normMag,binMaskref,binv,MaskPupil,rSize, ...
monitorSize,scrnIdx,coordType,abs_ang,MaxMask,plotMaskref,maskSelref);      

 %%% Register the reference:
 if measSimulated == 0
%     src.Exposure = tc0Exposure; % Change camera's exposure due to the 
                                                   % energy spreading: Experimental point for tc=0
    snap = getsnapshot(vid); % Real measurements
    wait(vid);
 else % measSimulated = 1
    snap = refMask; % "Simulated" measurements (the mask is saved)
    snap(snap==0) = 1; % Convert zeros to ones
    snap(isnan(snap)) = 0; % Converts nan's to zeros
 end
 
%%% Save the reference:
refImgPath = strcat(imgpartPath,'reference',infoDelim,'tc', ... 
                            num2str(tcref),infoDelim,'gl',num2str(glref));
refmeasfullpath =  strcat(refImgPath,dataformat);
% Explanation: imwrite(variables,directory+filename+extension)
imwrite(snap,refmeasfullpath); 

%%% Wait and close the figure:
pause(recordingDelay); % Displays the mask for "recordingDelay" seconds  
if measSimulated == 0
    wait(vid); % Waits until vid is not running or logging
end

if iscell(reffig)
    if ishandle(reffig{1})
      close(reffig{1}); % Closes the reference mask
    end
else
    if ishandle(reffig)
      close(reffig); % Closes the reference mask
    end
end
    

%% End of the measurements
src.Exposure = initialExposure; % Restore the initial camera's exposure
% Author: PhD student Jens de Pelsmaeker VUB B-PHOT 2018, Brussels, Belgium
% MATLAB built in:
t2_dt = datetime;
disp('Measurement finished'); disp(t2_dt)
time = t2_dt - t1_dt; % Relative difference between start and stop
disp('Measurement took:'); % datestr(time,'SS') ' seconds'])
disp(time);

%% End notification
N = 3; % Number of beeps for the measurement
f_EndBeeps(N,beepSound);

end