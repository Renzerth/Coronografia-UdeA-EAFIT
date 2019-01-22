%% General program parameters 
% Always executes internally where it is required

% Measurement directory
% meas_dir = '';
meas_dir = 'C:\Users\saple\Dropbox\DAVID-SAMUEL\2018-2\Work Life\1_Internship\2_B-Phot-VUB\1_SPR Biosensor\Measurements\October\just tests-no measuring\interrogation';
% cd(meas_dir); % Measurement directory

%% Material dependent function, change to lens or prism; Au or Ag
% open f_laser_intensity; % Adjust until normalization is correctly reached

%% In several scripts (to be constantly changed):
precision = 3; % Precision of displayed results: significative digits (3)

%%% System measurement parameters
polarization = 'TM'; % TE (90deg) or TM (180deg)
optical_component = 'prism4'; % prism/lens + number: component identifer
material = 'Au'; % Au gold Ag silver
medium = 'water'; % water, sol2%, sol5% (methanol)
pinhole = '40um'; % 100um % pinhole diameter. Default: 40um
dist = '42mm'; % distance of light source to component (indicated on the 
% horizontal syringe that is pasted from 40 to 60 mm). This distance
% shouldn't be changed but any modification may improve if the focus of the
% freeform lens is taken into account     

%%% Source parameter
source_current = 90; % Current in mA; 
                     % 100% laser intensity at 128mA (71.0373mW): 1.6183e8          

%%% Scan parameters
angle_or_intensity = 1; % angle (0) or intensity (1) interrogation

if angle_or_intensity == 0
    % Angle interrogation:
    startpos = 74; % Usually one only changes this
    angle_inc = 0.5; % add slight changes for Ag Au
    wait_inc = 2; % In seconds
    iterations = 20;
else
    % Intensity interrogation:
    startpos = 89; % Put the optimal point angle here, which is the 
                   % output of the "normalize" script. 
                   % Its value is the one of the variable 'max_ang'
    angle_inc = 0.1; % 0.1 or 0.2
    wait_inc = 0.7; % Waits less so that the measures are closer
    iterations = 10;
    startpos = startpos - angle_inc*iterations/2; % In order to scan at 
                                                  % both sides of the 
                                                  % optimal point
end

% Comments
% startpos: initial angle (check with f_readPos(s)); in degrees
% angle_inc: angle increment (degrees)
% wait_inc: should be 0.7 for 0.1 step % wait increment between 
%           measurements 2 for 0.5 degree step
% iterations: -LINUX: 100 takes 2m48 w 4 avg and 0.5 wait 
%             -WINDOWS: 3m19 (199s) w 4 averaging and 0.2 wait; ref:  20

endpos = startpos + iterations*angle_inc;
showpos = 1; % Every time the program is runned, show positions
if showpos == 1
disp(['Start position: ' num2str(startpos) 'deg; end position will ' ...
      'be: ' num2str(endpos) 'deg']);
end
% Example: startpos 80 with angle_inc 0.5 will arrive to 90 with 20 iterations
wait = 13; % 10 degrees is 12 sec; Wait while the rotor gets the laser on
           % the start position; As well one can cancel it ONLY on this 
           % time interval with ctrl+c

%% In several scripts (fixed):
level = 500; % Filtering level (bootstrap, empirically determined) 
             % 750 before, too high. Ref: 500. It deletes the Airy disks 
             % and only keeps the central part        
avg = 4; % Frames used for averaging: always 4 since it is enough and 
         % doesn't take too much time
         % The normalized reflectivity is found by summing all the pixels
         % on the image
           
%% Get_shot
avg_shot = 1; % One shot at a time
bootstrap = level; % Bootstrap is the level below which we clip the value 
                   % of a pixel  
show_parameters = 0; % No: 0; yes: 1

%% Intensity stability
% REMEMBER that startpos is used here: modify it before executing it
% Measurement:
lighting = 'lights_off'; % lights_off, lights_on
% Scan:
iterations_int_stab = 100; % 100 takes 1m39, 500 takes 8m17 (497 s)

%% Angle scanner
capt = 1; % Show captures (1) per angle
save_capt = 0; % Save the captures inside "shots" folder: for debugging
               % save_capt = 1 implies that capt is considered as 1.
               % (0) means deactivated
               
%% Normalization, general parameters (Normalization and reproducibility)
% These parameters apply for: normalize, normalize_meas and
% normalize_reprod

% Normalization calibration: norm_factor and f_laser_intensity
% Fix parameters FIRST on 'f_laser_intensity', norm_factor should ALWAYS
% equal one
norm_factor = 1; % Note: this factor helps more with "normalize_reprod" 
                 % than f_laser_intensity, but should NOT be changed
% open f_laser_intensity;

%%% Debugging: check how high is the plot when finding the correct
  %            normalization factor.
% xlim([0 21])
% ylim([0 2])     

keep_TM = 0; % Only saves TM and doesn't normalize with the TE
             % Just keep TM...NOT OK! Why not? Better then with TE
             % This is a "illegal" step. Doesn't do it (0 and morally 
             % correct); does it (1 and morally incorrect)
   
%%%  Normalization Method                              
% Change in case there are problems with the normalization
% Try to maintain normalization of norm_TM
% By diminishing influence of TE signal and placing around 1
norm_method = 1; % 1: normalize with TE, maxed at 1 (Default)
                 % 2: seems to give best performance. Influence of TE is
                 %    halved
                 % 3: do nothing
% The selected method used to be performed different for each program:
% normalize: method 1
% normalize_reprod: method 3
% normalize_meas: method 1


%% Normalize reproducibility (1 liquid): intensity interrogation
count = 2; % number of reproducibility measurements performed to normalize
           % here (TM/TE pairs)
           % This number MUST be bigger than 1 (minimum 2 measurements 
           % under same conditions and with in between "washing")
           
%% Reproducibility (1 liquid): intensity interrogation
word = 'repeatability'; % reproduciblity % repeatability
                        % Reproducibility: when I wash between measurements
                        % Repeatability: no wash between the measurements
descriptor= ' of 3 measurements over 3 minutes = '; % reproduciblity of N
                                                    % measurements over M 
                                                    % minutes
                                                     
%% Normalize measurement (until 6 liquids): intensity interrogation
% Manually put all solutions here and on the array (cell) of all the plots
sol1 = 'water';
sol2 = 'sol2%';
sol3 = 'sol5%';
medium_liquids = {sol1, sol2, sol3}; % Array of previous variables

%% Compare until 6 plots
Ksigma = 3; % Number of additional standard deviations that are multiplied
            % to the LOD so that it takes into account noise and give a 
            % more realistic estimation. It is a empirical constant that 
            % appears on the literature. 3 sigma, trust interval (99.7%)
XC = [ 0 1 2]; % Vector of different concentrations (percentages)
Kn = 1; % Coefficient of concetration to RI: [RIU]/[concentration units]
XRI = Kn*XC; % Refractive index as X vector    
residualBool = 0; % plot (1) or don't (0); Default: 0
errorBool = 1; % plot (1) or don't (0); Default: 1 
% 'rep_value' is taken from "reproducibility.m" and used in compare 6 or less plots
% In case it does not exist, define it here so that a reprod/repeat
% measruement is not needed to be made again (pre-taken value)
if ~exist('rep_value', 'var') % True if it does not exist
    rep_value = 0.0258; % This will only be applied if the script 
                        % "reproducibility.m" was not previously executed
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% OLD %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Refractive index difference = 5*10^-4; difference between each of the
% solutions, it is an approximation, the exact values ARE ON THE THESIS.
% It represents an almost-constant difference between H2O-sol2% and also
% between the sol2%-sol5% solutions
% It was stablished by trial/error as the lowest value that one could 
% assume to reliably differentiate/produce in the laboratory. It is the 
% standard error of making each solution and avoids measuring each solution
% on the Photonics Innovation Center.            
% ref_diff = 5e-4; % This constant is a good estimative used by Jens