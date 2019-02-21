% [imgfullpath] = f_AutomatMeasure(savetype,pathSep,dataformat,cameraPlane)

%% Automated measurements
showM = 0; % Don't show a fig in "PhaseMaskSel": this should always be 0.
totalImgs = ltcvect*lglvect; % Number of images to be taken
expImgs = cell(1,totalImgs); % Cell with the experimental images
MeasInfo = expImgs; % Same initialization as expImgs
idxgral = 1; % Initialization of the general index that runs 
             % on: [1,totalImgs]
% fileFormat = '.bmp'; % OLD: save each image

%% Measurements initialization
% Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
pause(wait) % Seconds before measuring as a safety measurement
t1_dt = datetime; % store time
disp('Measurement started'); disp(t1_dt)
tit = 'Displayed phase mask';
imgpath = strcat(DatalogDir,pathSep,cameraPlane,'_'); % More information 
% will be concatenated for a full path of the measured images inside the
% next "for" loops

%% Measurements
for idxtc = 1:ltcvect 
  for idxgl = 1:lglvect
      
    %% Generate the phase mask
    tc = tcvect(idxtc); % Specific tc for this iteration
    gl = glvect(idxgl); % Specific gl for this iteration
    PhaseMaskSel; % Selects a phase mask to display
    wrappedMask = f_mask_circ_angle_gl(r,mask,binMask,glphi,mingl, ...
                                       maxgl,levShft);
                                   
    %% Display the phase mask on the SLM
    slmhfig = f_fig_maskSLM(x,y,r,mask,gl,glphi,mingl,maxgl,levShft, ... 
                            abs_ang,binMask,monitorSize,plotMask);
   
    %% Record a snapshot
    if measSimulated == 0
        snap = getsnapshot(vid); % Real measurements
    else % measSimulated = 1
        snap = wrappedMask; % "Simulated" measurements (the mask is saved)
    end
    expImgs{idxgral} = snap; % An extructure for future use (?)
    A = expImgs{idxgral}; % A variable for the save function

    %% Saving the measurement
    tcstr = strcat('tc_',num2str(tcvect(idxtc))); 
    glstr = strcat('gl_' num2str(glvect(idxgl)));
    MeasInfo{idxgral} = [tcstr '_' glstr]; % Dataname for each experimental
                                           % data
    if savetype == 1  % .mat format
        imgfullpath = [imgpath MeasInfo{idxgral}];
        % save(directory+filename,variables) % ,'-append'
        save(imgfullpath,'A'); % Save as .matd
    else % savetype = 2. dataformat is used
        imgfullpath = [imgpath MeasInfo{idxgral} dataformat];
        % imwrite(variables,directory+filename+extension)
        imwrite(expImgs{idxgral}, imgfullpath); 
    end
        
    %% Displaying the measurement
% Author: PhD student Jens de Pelsmaeker VUB B-PHOT 2018, Brussels, Belgium
    % The numbers after 'position' were empirically obtained
    fig = figure('units','normalized','position',[1/10 1/10 1/3 1/2]);
    imagesc(snap); % normalized
    colorbar; title(strcat('Camera image: ',filename));

    % The numbers after 'position' were empirically obtained
    showmask = 1;
    pcfig = f_fig_maskPCscreen(x, y, wrappedMask, tit, gl, showmask);
    set(pcfig,'units','normalized','position',[5/10 1/10 1/3 1/2]);
     
    %% Preparation for a new measurement iteration          
    stridxgral = num2str(idxgral); strtotalImgs = num2str(totalImgs);
    disp(strcat(stridxgral,' out of ',strtotalImgs,' images recorded'));
    idxgral = idxgral + 1; % The general index increases   
    pause(recordingDelay); % Displays the mask for "recordingDelay" seconds   
                           % This time is also important so that the camera
                           % bus doesn't overload 
    close(pcfig); close(fig); close(slmhfig); % Close the displayed figures
  end
end
% MATLAB 2018b: disp(newline); MATLAB 2016: disp(char(10)) 
% Aeasthetic reasons: Not needed in MATLAB 2016

%% End of the measurements
% Author: PhD student Jens de Pelsmaeker VUB B-PHOT 2018, Brussels, Belgium
% MATLAB built in:
t2_dt = datetime;
disp('Measurement finished'); disp(t2_dt)
time = t2_dt - t1_dt;
disp('Measurement took:'); % datestr(time,'SS') ' seconds'])
disp(time)
% end