%% Automated measurements
showM = 0; % Don't show a fig in "PhaseMaskSel": this should always be 0.
totalImgs = ltcvect*lglvect; % Number of images to be taken
expImgs = cell(1,totalImgs); % Cell with the experimental images
MeasInfo = expImgs; % Same initialization as expImgs
idxgral = 1; % General index that will be on [1,totalImgs]
% fileFormat = '.bmp'; % OLD: save each image

for idxtc = 1:ltcvect 
  for idxgl = 1:lglvect
    tc = tcvect(idxtc); % Specific tc for this iteration
    gl = glvect(idxgl); % Specific gl for this iteration
    cd(analysDir); % Moves to the Analysis directory
    PhaseMaskSel; % Selects a phase mask to display
    f_fig_maskSLM(x,y,r,mask,m,n,a,b,gl,abs_ang,binMask,plotMask);
    pause(recordingDelay); % Displays the mask for "recordingDelay" seconds   
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
    % imwrite(expImgs{idxgral},[MeasInfo{idxgral} fileFormat]); % OLD: save each image % Saves the last shown figure
    % Here, the masks are saved, but later the images should
    % be saved as an input of the algorithm to be processed
    % f_CameraShot();
    idxgral = idxgral + 1; % The general index increases   
  end
end
cd(analysDir); % Go to script Directory