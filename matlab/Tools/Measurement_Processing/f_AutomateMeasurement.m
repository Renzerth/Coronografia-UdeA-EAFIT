function [] = f_AutomateMeasurement(Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,...
phiPC,phaseValues,mingl,maxgl,levShft,s,ph0,p,WsizeRatio,L,f_FR,bcst,period,T0, ...
frkTyp,Aalpha,Angalp,Angbet,z_coeff,a,frac,pupil,ZernikeSize,disp_wrap,...
plot_z, normMag,binMask,binv,MaskPupil,monitorSize,scrnIdx,coordType,abs_ang, ...
MaxMask,maskSel,ltcvect,lglvect,wait,DatalogDir,dataformat,pathSep,cameraPlane,...
tcvect,glvect,measSimulated,recordingDelay) 
% Plots phase masks on the Fourier plane of the vortex coronagraph and
% takes images of either its Lyot or PSF plane

%% Automated measurements
totalImgs = ltcvect*lglvect; % Number of images to be taken
expImgs = cell(1,totalImgs); % Cell with the experimental images
MeasInfo = cell(1,totalImgs); % Same initialization as expImgs
idxgral = 1; % Initialization of the general index that runs 
             % on the range: [1,totalImgs]

%% Measurements initialization
% Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
pause(wait) % Seconds before measuring as a safety measurement
t1_dt = datetime; % store time
disp('Measurement started'); disp(t1_dt)
% tit = 'Displayed phase mask'; Not used right now: title for the "Display
% the mask on the PC" section
imgPath = strcat(DatalogDir,pathSep,cameraPlane,'_'); % More information 
% will be concatenated for a full path of the measured images inside the
% next two "for" loops

%% Measurements
for idxtc = 1:ltcvect 
  for idxgl = 1:lglvect
      
    %% Generate the phase mask and display it on the SLM
    tc = tcvect(idxtc); % Specific tc for this iteration
    gl = glvect(idxgl); % Specific gl for this iteration
    plotMask = 2; % Select SLM
    [X,Y,r,phi] = f_SelectCoordinates(Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,...
                                      phiPC,plotMask);
    [~,wrapMask,slmfig] = f_PlotSelectedMask(X,Y, ...
r,phi,gl,phaseValues,mingl,maxgl,levShft,tc,s,ph0,p,WsizeRatio,L,f_FR,bcst,period,T0, ...
frkTyp,Aalpha,Angalp,Angbet,z_coeff,a,frac,pupil,ZernikeSize,disp_wrap,...
plot_z, normMag,binMask,binv,MaskPupil,monitorSize,scrnIdx,coordType,abs_ang, ...
MaxMask,plotMask,maskSel);  
  
    %% Record a snapshot
    if measSimulated == 0
        snap = getsnapshot(vid); % Real measurements
    else % measSimulated = 1
        snap = wrapMask; % "Simulated" measurements (the mask is saved)
    end
    expImgs{idxgral} = snap; % An extructure for future use (?)
    % A = expImgs{idxgral}; % A variable for the save function

    %% Saving the measurement
    tcstr = strcat('tc_',num2str(tcvect(idxtc))); 
    glstr = strcat('gl_',num2str(glvect(idxgl)));
    MeasInfo{idxgral} = [tcstr '_' glstr]; % Dataname for each experimental
                                           % data
    imgfullpath = strcat(imgPath,MeasInfo{idxgral},dataformat);
    % imwrite(variables,directory+filename+extension)
    imwrite(expImgs{idxgral}, imgfullpath); 
        
    %% Displaying the camera on the PC
% Author: PhD student Jens de Pelsmaeker VUB B-PHOT 2018, Brussels, Belgium
    % The numbers after 'position' were empirically obtained
    camfig = figure('units','normalized','position',[1/10 2/10 1/3 1/2]);
    imagesc(snap); % normalized
    colorbar; title(strcat('Camera image: ',cameraPlane));
    % The numbers after 'position' were empirically obtained
    
    %% Display the mask on the PC
    plotMask = 1; % Select PC
    [X,Y,r,phi] = f_SelectCoordinates(Xslm,Yslm,rSLM,phiSLM,Xpc,Ypc,rPC,...
                                      phiPC,plotMask);                             
    [~,~,pcfig] = f_PlotSelectedMask(X,Y, ...
r,phi,gl,phaseValues,mingl,maxgl,levShft,tc,s,ph0,p,WsizeRatio,L,f_FR,bcst,period,T0, ...
frkTyp,Aalpha,Angalp,Angbet,z_coeff,a,frac,pupil,ZernikeSize,disp_wrap,...
plot_z, normMag,binMask,binv,MaskPupil,monitorSize,scrnIdx,coordType,abs_ang, ...
MaxMask,plotMask,maskSel);  
    set(pcfig,'units','normalized','position',[5/10 2/10 1/3 1/2]);
     
    %% Preparation for a new measurement iteration          
    stridxgral = num2str(idxgral); strtotalImgs = num2str(totalImgs);
    disp(strcat(stridxgral,' out of ',{' '},strtotalImgs, ...
         ' images recorded'));
    idxgral = idxgral + 1; % The general index increases   
    pause(recordingDelay); % Displays the mask for "recordingDelay" seconds   
                           % This time is also important so that the camera
                           % bus doesn't overload 
    close(pcfig); close(camfig); close(slmfig); % Close the displayed figures
  end
end
% MATLAB 2018b: disp(newline); MATLAB 2016: disp(char(10)) 
% Aeasthetic reasons: Not needed in MATLAB 2016

%% Store all measurements in a .mat file
imgfullpath = strcat(imgPath,'allmeas'); % Saves the cell
% save(directory+filename,variables) % ,'-append'
save(imgfullpath,'expImgs'); % Save as .mat

%% End of the measurements
% Author: PhD student Jens de Pelsmaeker VUB B-PHOT 2018, Brussels, Belgium
% MATLAB built in:
t2_dt = datetime;
disp('Measurement finished'); disp(t2_dt)
time = t2_dt - t1_dt; % Relative difference between start and stop
disp('Measurement took:'); % datestr(time,'SS') ' seconds'])
disp(time)
% end