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
disp('start '); disp(t1_dt)

%% Measurements
for idxtc = 1:ltcvect 
  for idxgl = 1:lglvect
    tc = tcvect(idxtc); % Specific tc for this iteration
    gl = glvect(idxgl); % Specific gl for this iteration
    cd(analysDir); % Moves to the Analysis directory
    PhaseMaskSel; % Selects a phase mask to display
    
    
    %% Show frame
    % Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
%     figure('Position',[250 100 700 500] ); imagesc(SingleFrame); 
%     colorbar; title(['Camera image: ' filename]);
%     figure('Position',[1000 100 700 500] ); imagesc(log(SingleFrame)); 
%     colorbar; title(['LOG camera image: ' filename])
    
    % OR:
    f_fig_maskSLM(x,y,r,mask,m,n,a,b,gl,abs_ang,binMask,plotMask);
    
    
    
    pause(recordingDelay); % Displays the mask for "recordingDelay" seconds   
                           % This time is also important so that the camera
                           % bus doesn't overload
    tcstr = ['tc_' num2str(tcvect(idxtc))]; 
    glstr = ['_gl_' num2str(glvect(idxgl))];
    MeasInfo{idxgral} = [tcstr glstr]; % Dataname for each experimental data
    cd(DatalogDir); % Goes to the data log directory (specific measurement
                    % folder)
    wrappedMask = f_circularPupil_maskAngle(r,mask,binMask); 
    snap = wrappedMask; % "Simulated" measurements (the mask is saved)
    % snap = getsnapshot(vid); % Real measurements
    expImgs{idxgral} = snap;
    A = expImgs{idxgral};
    save(MeasInfo{idxgral},'A'); % save(filename,variables,'-append')
    
    % OLD: save each image % Saves the last shown figure
    % imwrite(expImgs{idxgral},[MeasInfo{idxgral} fileFormat]); 
    % Here, the masks are saved, but later the images should
    % be saved as an input of the algorithm to be processed
    % f_CameraShot();
    idxgral = idxgral + 1; % The general index increases   
  end
end
cd(analysDir); % Go to script Directory

%% End of the measurements
% Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
set(0,'DefaultFigureVisible','on'); 
% MATLAB built in:
t2_dt = datetime;
disp('stop '); disp(t2_dt)
time = t2_dt - t1_dt;
disp('It took '); % datestr(time,'SS') ' seconds'])
disp(time)