function [refmeasfullpath] = f_AutomateMeasurement(vid,Xslm,Yslm,rSLM, ...
phiSLM,Xpc,Ypc,rPC,phiPC,s,ph0,p,WsizeRatio,L,f_FR,bcst,period,T0, ...
frkTyp,Aalpha,Angalp,Angbet,z_coeff,z_a,z_pupil,z_disp_wrap,z_plot, ...
normMag,binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx,coordType, ...
abs_ang,MaxMask,maskSel,ltcvect,lglvect,totalImgs,wait,imgpartPath, ...
dataformat,imgformat,measfullpath,infoDelim,cameraPlane,tcvect,glvect, ...
measSimulated,recordingDelay)
% Plots phase masks on the Fourier plane of the vortex coronagraph and
% takes images of either its Lyot or PSF plane
%
% Inputs:
%  vid: video input object
%
% Ouputs:

%% Automated measurements
expImgs = cell(1,totalImgs); % Cell with the experimental images
MeasInfo = cell(1,totalImgs); % Same initialization as expImgs
idxgral = 1; % Initialization of the general index that runs 
             % on the range: [1,totalImgs]

%% Measurements initialization
% Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
pause(wait); % Seconds before measuring as a safety measurement
t1_dt = datetime; % store time
disp('Measurement started'); disp(t1_dt)

%% Measurements
for idxtc = 1:ltcvect 
  for idxgl = 1:lglvect
      
   %% Generate the phase mask and display it on the SLM
    tc = tcvect(idxtc); % Specific tc for this iteration
    phaseValues = glvect(idxgl); % Specific phase values for this iteration
    phaseValues = linspace(0,2*pi,phaseValues);
    plotMask = 2; % Select SLM for plotting
    [X,Y,r,phi] = f_SelectCoordinates(Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,...
                                      phiPC,plotMask);
    [~,wrapMaskslm,slmfig,~] = f_PlotSelectedMask(X,Y,r,phi,phaseValues,tc, ...
    s,ph0,p,WsizeRatio,L,f_FR,bcst,period,T0,frkTyp,Aalpha,Angalp, ...
    Angbet,z_coeff,z_a,z_pupil,z_disp_wrap,z_plot,normMag,binMask,binv, ...
    MaskPupil,rSize,monitorSize,scrnIdx,coordType,abs_ang,MaxMask, ...
    plotMask,maskSel);      
  
   %% Record a snapshot
    if measSimulated == 0
        snap = getsnapshot(vid); % Real measurements
    else % measSimulated = 1
        snap = wrapMaskslm; % "Simulated" measurements (the mask is saved)
    end
   
   %% Displaying the camera on the PC
    % The numbers after 'position' were empirically obtained
    camfig = figure('units','normalized','position',[1/12 2/10 3/7 1/2]);
                                               % [left bottom width height]
    imagesc(snap); % normalized
    colorbar; title(strcat('Camera image:',{' '},cameraPlane));
    pax = gca; pax.FontSize = 16; % Font size
    % The numbers after 'position' were empirically obtained
    
    %% Display the mask on the PC
    plotMask = 1; % Select PC
    [X,Y,r,phi] = f_SelectCoordinates(Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,...
                                      phiPC,plotMask);                             
    [~,wrapMaskpc,pcfig,~] = f_PlotSelectedMask(X,Y,r,phi,phaseValues, ...  % TEMPORARLY NOT BEING USED
    tc,s,ph0,p,WsizeRatio,L,f_FR,bcst,period,T0,frkTyp,Aalpha,Angalp, ...
    Angbet,z_coeff,z_a,z_pupil,z_disp_wrap,z_plot,normMag,binMask,binv, ...
    MaskPupil,rSize,monitorSize,scrnIdx,coordType,abs_ang,MaxMask, ...
    plotMask,maskSel);     
    set(pcfig,'units','normalized','position',[6/11 2/10 3/7 1/2]);

    
   %% Saving the measurement
    if measSimulated == 0
         expImgs{idxgral} = snap; % An extructure with all the images
    else % measSimulated = 1
         % If you want to simulate with the shifted mask, put wrapMaskslm
         expImgs{idxgral} = wrapMaskslm; % Saves the mask
    end
    
    tcstr = strcat('tc',infoDelim,num2str(tcvect(idxtc))); 
    glstr = strcat('gl',infoDelim,num2str(glvect(idxgl)));
    MeasInfo{idxgral} = strcat(tcstr,infoDelim,glstr); % Dataname for each 
                                           % experimental data
    imgfullpath = strcat(imgpartPath,MeasInfo{idxgral});
    if measSimulated == 0
        % Explanation: imwrite(variables,directory+filename+extension)
        imwrite(expImgs{idxgral}, strcat(imgfullpath,dataformat)); 
    else
        % Explanation: saveas(variable,directory+filename,extension)
        saveas(gcf,imgfullpath,imgformat); % Saves the last shown figure
        % Other options: savefig and print.
    end
     
    %% Preparation for a new measurement iteration          
    stridxgral = num2str(idxgral); strtotalImgs = num2str(totalImgs);
    disp(strcat(stridxgral,' out of ',{' '},strtotalImgs, ...
         ' images recorded'));
    idxgral = idxgral + 1; % The general index increases   
    pause(recordingDelay); % Displays the mask for "recordingDelay" seconds   
                           % This time is also important so that the camera
                           % bus doesn't overload 
    close(pcfig); close(camfig); close(slmfig); % Close the shown figures
  end
end

%% Store all measurements in a .mat file
% Explanation: save(directory+filename,variables) % ,'-append'
save(measfullpath,'expImgs','MeasInfo'); % Save as .mat

%% Reference measurement
% Without tc and for 256 gray levels (white background)
% This is the non-coronagraphic PSF reference: the system response without
% phase masks

%%% Generate reference mask:
tcref = 0;
phaseValuesref = [0 255]; % [0 255]: white; [0 1]: black
phaseValuesref = phaseValuesref*2*pi/255; % Conversion from gray-levels to 
                                 % phase levels      
pref = 0;
binMaskref = 1;
maskSelref = 1; % LG mask
plotMaskref = 2; % Select SLM for plotting
WsizeRatioref = 100;
[Xref,Yref,rref,phiref] = f_SelectCoordinates(Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,...
                                  phiPC,plotMaskref);
[~,wrapMaskslmref,~] = f_PlotSelectedMask(Xref,Yref,rref,phiref,phaseValuesref,tcref, ...
s,ph0,pref,WsizeRatioref,L,f_FR,bcst,period,T0,frkTyp,Aalpha,Angalp, ...
Angbet,z_coeff,z_a,z_pupil,z_disp_wrap,z_plot,normMag,binMaskref,binv, ...
MaskPupil,rSize,monitorSize,scrnIdx,coordType,abs_ang,MaxMask, ...
plotMaskref,maskSelref);      

%%% Register the reference:
 if measSimulated == 0
    snap = getsnapshot(vid); % Real measurements
else % measSimulated = 1
    snap = wrapMaskslmref; % "Simulated" measurements (the mask is saved)
 end
close(gcf); % Closes the reference mask
 
% Save the reference:
refImgPath = strcat(imgpartPath,'reference',infoDelim,'tc',num2str(tcref),infoDelim,'gl,',num2str(256));
refmeasfullpath =  strcat(refImgPath,dataformat);
% Explanation: imwrite(variables,directory+filename+extension)
imwrite(snap,refmeasfullpath); 

%% End of the measurements
% Author: PhD student Jens de Pelsmaeker VUB B-PHOT 2018, Brussels, Belgium
% MATLAB built in:
t2_dt = datetime;
disp('Measurement finished'); disp(t2_dt)
time = t2_dt - t1_dt; % Relative difference between start and stop
disp('Measurement took:'); % datestr(time,'SS') ' seconds'])
disp(time)
% end