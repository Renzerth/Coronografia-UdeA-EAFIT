%% Principal program: SPR measurements
% Jens (major work) and Samuel (comments and organization, Jul-Dec 2018); 
% Tatevik: ideas on calculating the LOD and generalities
%
% Information:
%  Check "SPR_set-up procedures.docx" for further info
%  Do NOT press RUN, better execute which step when ready
%  The 'open' command is used for scripts that must be always edited
%  Ctrl + Enter to execute a part between the sections splitted by '' %% ''
%  If the camera doesn't work: connect and disconnect it and reboot linux
%  set(0,'DefaultFigureVisible','off'); 'on' or 'off' for all plots
%  Motor off emergecy? Disconnect DC alimentation
%  Polarizer positions in degrees:
%   180: TM (p-polarized, plasmon couples with it)
%   90: TE (should be constant, useful for normalization)
% In Linux, the degree symbol is not identified so it will appear as a
% question mark

% Missing:
% If the camera works on windows, then it can be used both on Windows and
% Linux

% Questions:
% open f_power_normalizer: I do not know what does it do.
% compare_3_plots: lines 119 and 128

%% 0_Always execute
clc; clear; close all; % Initialization
% CHECK CURRENT DIRECTORY AND PROGRAM VERSION:
windows_or_linux = 1; % Windows(0) or Linux(1)
if windows_or_linux == 1
    cd '/home/jens/Documents/Jens-Samuel_SPR-Matlab/28'; % Linux
else
    cd 'C:\Users\saple\Dropbox\DAVID-SAMUEL\2018-2\Work Life\1_Internship\2_B-Phot-VUB\1_SPR Biosensor\Matlab-Jens_Samuel-SPR setup control\28'; % Samuel's mac
end
% Add directories
addpath('Functions'); % f_add_paths is inside this folder
f_add_paths; % Adds all folders inside Jens-Samuel_SPR/Nth_version
% In order to restore back default paths, type: restoredefaultpath
scripts_dir = pwd; % Default script's directory: save it when in there
                   % (not used so far)

%% 1_Initial-settings
% In windows, the COMx needs to be changed inside this script depending on
% the used port
Initiate_settings; % Turns camera on and all the needed variables
% 'Get_camera_on' is used on this script
% Remember to leave the preview open
% 
% Usually these warnings shows up:
% Warning 1: The LineSource property could not be set.  Check the current 
%            value before continuing your work. 
% Solved by leaving the preview open at all the times, otherwise the
% measurements become slow. In order to open it again, execute 
% "Get_camera_on" on the Command Window so that preview(vid) is executed
%
% Warning 2: MATLAB has disabled some advanced graphics rendering features 
%            by switching to software OpenGL.
% Will delay the feedback of the camera only.
% 
% Warning 3: No Image Acquisition adaptors found. Image acquisition 
%            adaptors may be available as downloadable support packages. 
%            Open Support. Package Installer to install additional vendors.
% Solved by  making a reboot of windows

%% General program parameters 
Parameters; % Always executes internally where it is required
open Parameters; % Keep open always

%% 2_Tests
% When in doubts, when aligning, use it!
Get_shot; % Takes a photo, shows a figure of it and saves it as shot.png

% REMEMBER that startpos is used here: modify it before executing it
intensity_stability; % "power_normalizer" is used inside it.
                        % Checks laser source stability, should flucutate
                        % less than 1% is ok; Takes about 2min
                        % Check it at least 3 times before continuing
                        % Results given:
                        %  1)Each measurement: punctual standard deviation;
                        %    it is the std of the average of 3 images,
                        %    represents the CCD stability when taking images
                        %  2)Source power: represents the laser stability
                        %    standard deviation
                        % Normal warning: Warning: Directory already exists. 
                        % Normally takes 2 minutes                        
                        
f_movePos(s,-2); pause(3); f_movePos(s,2);  %  test that the rotator works 
                                            % with a back and fort!: pause of 2s

%% 3_Measurement: used for ALL types of measurements
%%% Angle interrogation
f_setPos(s,startpos); pause(2); Angle_scanner; % Fast switching between
                                               % TE/TM; 2s pause

%%% Intensity interrogation                                                  
f_setPos(s,startpos);                             
Angle_scanner; % Scans at different angles
               % Inside: f_laser_intensity, f_LD_c2p and 
               % f_power_normalizer

% The reflectivity is found by first finding the sum of all the pixel 
% values (energy of the acquired image with the bootstrap already applied)
% of a specific snapshot for which one obtains a vector with the energy of 
% each snapshot repetition (usually 4 consecutive photos are taken in a row 
% per measurement for a specific angle). Finally, the mean and the standard 
% deviation are used and this procedure is repeated for each angle so that
% one obtains the 1D plot of reflectivity vs. angle. The average of the 4 
% images corresponds at the snapshot that is shown every time one previews 
% the camera with get_snapshot.
 
                  
%%%%%%%%% 4_Normalization   
%% Angle interrogation: independent measurement (finds the optimal point)
normalize; % 10 degrees Au or 5 degrees Ag
           % Just to check one solution (liquid) measurement of TE/TM and 
           % shows it with the noise processed (minimized wrt the TE
           % measurement, that doesn't present SPR)
           % f_LD_c2p, f_laser_intensity are used
           % Note: an error in "f_laser_intensity" wouldn't be shown but 
           % won't let "normalize" execute
           % Note: this script is similar to normalize_reprod and
           % normalize_meas. This has old pieces of code commented in case
           % they are needed
open f_laser_intensity; % Adjust until normalization is correctly reached
           
%% Intensity interrogation: reproducibility or repeatability:      
% 1 solution, 3 measurements are enough, 1 or 2 degrees scanned
normalize_reprod; % Checks reproducibility for one solution (liquid)
                  % ''f_LD_c2p'', ''f_laser_intensity'' are used
                  
% Reproducibility: when I wash between the solutions
% Repeatability: no wash between the solutions

% The ''reproducibility'' program needs: ''normalize_reprod''
reproducibility; % Checks reproducibility for one solution (liquid)
                 % Also finds the optimal point for reprod. and repeatab.

%% Measurement of SPR: Intensity interrogation 
% Solutions: two to six. 
% 3 solutions are enough, 1 or 2 degrees scanned around the optimal point
normalize_meas; % Measurement of 3 different solutions

%% 5_Comparison (directly related with: normalize_meas; measurement of SPR)
% The ''compare_until_6_plots'' program needs: ''normalize_repeatab''
% ''reproducibility'' and ''normalize_meas''; to had been previously done
% Solutions: two to six
compare_until_6_plots; % Includes Limit Of Detection (LOD)

%% Functions
f_movePos(s,-1); % Relative movement, degrees
f_setPos(s,95); % Absolute movement, degrees. Range: [10,105] TAKE CARE
                % Camera's field of view: [80,91] ~ 10 degrees
f_readPos(s); % Shows current absolute position

f_laser_intensity; % Compensates normalization depending on laser power
                   % and other system parameters

f_LD_c2p(source_current); % Transforms current to power with a previously
                          % measured data (performs a linear 
                          % inter/extra-polation)
f_power_normalizer; % Used in angle scanner

%% Not used
compare_2_plots; % Not used; doesn't calculate LOD
compare_5_plots; % Not used; doesn't calculate LOD
calibrate_noise; % Not used; checks the power used, without a source. 
                   % This is the dark current, it used to be implemented on 
                   % f_power_normalize as the noise_per_pixel variable for
                   % the background removal feature, now it's just bootstrap

%% 6_Termination
f_setPos(s,95); % Default position 
Terminate_settings; % Clears variables and closes all; deactivates camera 
                    % and laser rotation mechanism               

% Type exit on the Linux terminal                    