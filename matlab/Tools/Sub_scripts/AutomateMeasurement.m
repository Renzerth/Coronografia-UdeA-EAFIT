% [imgfullpath] = f_AutomateMeasurement(savetype,pathSep,dataformat,cameraPlane)
% Plots phase masks on the Fourier plane of the vortex coronagraph and
% takes images of either its Lyot or PSF plane
%% Automated measurements
showM = 0; % Don't show a fig in "PhaseMaskSel": this should always be 0.
totalImgs = ltcvect*lglvect; % Number of images to be taken
expImgs = cell(1,totalImgs); % Cell with the experimental images
MeasInfo = cell(1,totalImgs); % Same initialization as expImgs
idxgral = 1; % Initialization of the general index that runs 
             % on: [1,totalImgs]
% fileFormat = '.bmp'; % OLD: save each image

%% Measurements initialization
% Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
pause(wait) % Seconds before measuring as a safety measurement
t1_dt = datetime; % store time
disp('Measurement started'); disp(t1_dt)
tit = 'Displayed phase mask';
imgPath = strcat(DatalogDir,pathSep,cameraPlane,'_'); % More information 
% will be concatenated for a full path of the measured images inside the
% next "for" loops
plotMask = 2; % Always plot on the SLM



%% Measurements
for idxtc = 1:ltcvect 
  for idxgl = 1:lglvect
      
    %% Generate the phase mask and display it on the SLM
    tc = tcvect(idxtc); % Specific tc for this iteration
    gl = glvect(idxgl); % Specific gl for this iteration
    [mask,maskName] = f_PlotSelectedMask(X,Y,r,phi,gl,glphi,mingl, ...
    maxgl,levShft,tc,s,ph0,p,W,L,f_FR,bcst,period,T0,frkTyp,Aalpha, ...
    Angalp,Angbet,z_coeff,a,frac,pupil,sSize,disp_wrap,plot_z,normMag, ...
    binMask,binv,monitorSize,scrnIdx,coordType,abs_ang,plotMask,maskSel);  
     
    %% Record a snapshot
    if measSimulated == 0
        snap = getsnapshot(vid); % Real measurements
    else % measSimulated = 1
        snap = wrappedMask; % "Simulated" measurements (the mask is saved)
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
        
    %% Displaying the measurement
% Author: PhD student Jens de Pelsmaeker VUB B-PHOT 2018, Brussels, Belgium
    % The numbers after 'position' were empirically obtained
    fig = figure('units','normalized','position',[1/10 1/10 1/3 1/2]);
    imagesc(snap); % normalized
    colorbar; title(strcat('Camera image: ',filename));

    % The numbers after 'position' were empirically obtained
    showmask = 1;
    pcfig = f_ProjectMaskPC(x, y, wrappedMask, tit, gl, showmask);
    set(pcfig,'units','normalized','position',[5/10 1/10 1/3 1/2]);
     
    %% Preparation for a new measurement iteration          
    stridxgral = num2str(idxgral); strtotalImgs = num2str(totalImgs);
    disp(strcat(stridxgral,' out of ',{' '},strtotalImgs, ...
         ' images recorded'));
    idxgral = idxgral + 1; % The general index increases   
    pause(recordingDelay); % Displays the mask for "recordingDelay" seconds   
                           % This time is also important so that the camera
                           % bus doesn't overload 
    close(pcfig); close(fig); close(slmhfig); % Close the displayed figures
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