% [imgfullpath] = f_AutomatMeasure(savetype,pathSep,dataformat)

%% Automated measurements
showM = 0; % Don't show a fig in "PhaseMaskSel": this should always be 0.
totalImgs = ltcvect*lglvect; % Number of images to be taken
expImgs = cell(1,totalImgs); % Cell with the experimental images
MeasInfo = expImgs; % Same initialization as expImgs
idxgral = 1; % Initialization of the general index that runs 
             % on: [1,totalImgs]
% fileFormat = '.bmp'; % OLD: save each image

%% Measurements initializing
% Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
pause(wait) % Seconds before measuring as a safety measurement
t1_dt = datetime; % store time
disp('Measurement started'); disp(t1_dt)

%% Measurements
for idxtc = 1:ltcvect 
  for idxgl = 1:lglvect
    tc = tcvect(idxtc); % Specific tc for this iteration
    gl = glvect(idxgl); % Specific gl for this iteration
    PhaseMaskSel; % Selects a phase mask to display
    
    
    %% Show frame
% Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
%     figure('Position',[250 100 700 500] ); imagesc(SingleFrame); 
%     colorbar; title(['Camera image: ' filename]);
%     figure('Position',[1000 100 700 500] ); imagesc(log(SingleFrame)); 
%     colorbar; title(['LOG camera image: ' filename])
    
    % Display on the SLM
    f_fig_maskSLM(x,y,r,mask,gl,glphi,mingl,maxgl,levShft,abs_ang, ... 
                  binMask,monitorSize,plotMask);

    
    pause(recordingDelay); % Displays the mask for "recordingDelay" seconds   
                           % This time is also important so that the camera
                           % bus doesn't overload
    tcstr = ['tc_' num2str(tcvect(idxtc))]; 
    glstr = ['_gl_' num2str(glvect(idxgl))];
    MeasInfo{idxgral} = [tcstr '_' glstr]; % Dataname for each experimental
                                           % data
    wrappedMask = f_mask_circ_angle_gl(r,mask,binMask,glphi,mingl, ...
                                       maxgl,levShft);
    if measSimulated == 0
        % snap = getsnapshot(vid); % Real measurements
    else % measSimulated = 1
        snap = wrappedMask; % "Simulated" measurements (the mask is saved)
    end
    expImgs{idxgral} = snap;
    A = expImgs{idxgral};
    if savetype == 1  % .mat format
        imgfullpath = [DatalogDir pathSep MeasInfo{idxgral}];
        % save(directory+filename,variables) % ,'-append'
        save(imgfullpath,'A'); % Save as .mat
    else % savetype = 2. dataformat is used
        imgfullpath = [DatalogDir pathSep MeasInfo{idxgral} dataformat];
        % imwrite(variables,directory+filename+extension)
        imwrite(expImgs{idxgral}, imgfullpath); 
    end
    % OLD: save each image % Saves the last shown figure
    % imwrite(expImgs{idxgral},[MeasInfo{idxgral} fileFormat]); 
    % Here, the masks are saved, but later the images should
    % be saved as an input of the algorithm to be processed
    % f_CameraShot();
    stridxgral = num2str(idxgral); strtotalImgs = num2str(totalImgs);
    disp([stridxgral ' out of ' strtotalImgs ' images recorded']);
    idxgral = idxgral + 1; % The general index increases   
  end
end
disp(newline); % Aeasthetic reasons

%% End of the measurements
% Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
% MATLAB built in:
t2_dt = datetime;
disp('Measurement finished'); disp(t2_dt)
time = t2_dt - t1_dt;
disp('Measurement took: '); % datestr(time,'SS') ' seconds'])
disp(time)
% end