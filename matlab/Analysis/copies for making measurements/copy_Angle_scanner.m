% Used for intensity and angle interrogation, so a generic measurement on a
% specified range is done

close all
clearvars -except s vid vid_device snap_x snap_y max_ang; % clc
Parameters; % Always executes internally where it is required

%% Initialize folder
cd(meas_dir); % Measurement directory: always saves there in order to avoid 
              % confusions

%% Shot capturing at each angle
set(0,'DefaultFigureVisible','on');  % Stop displaying figures: 'off' or 'on'

% Used parameters:
% wait
% capt
% % Measurement parameters
% level
% optical_component
% material
% polarization
% medium
% source_current
% pinhole
% dist
% 
% Scan parameters:
% avg
% angle_inc
% wait_inc
% iterations
% startpos

% Choose window + calibration
if ~exist('snap_y','var'); snap_y = ':'; end
if ~exist('snap_x','var'); snap_x = ':'; end

noise = f_power_normalizer(snap_x,snap_y); % It is saved but not used here

%% Background operations before start
endpos = startpos + angle_inc*iterations;
add = [polarization '_' num2str(startpos) '-' num2str(endpos) '_' dist '_' optical_component '_' material '_' medium '_' pinhole '_' num2str(source_current) 'mA' ] ;
source_power = f_LD_c2p(source_current);

%% Initialize folder
cd(meas_dir); % Measurement directory
fldr_name = datestr(now);
fldr_name(fldr_name==':') = 'u'; 
fldr_name = fldr_name(1:end-3);
fldr_name  = [add ' (' fldr_name ')'];
mkdir(fldr_name); cd(fldr_name); 
if save_capt
 mkdir shots; % To store the shots taken
end

%% Initialize variables
counter = 0;
save('counter','counter')
angle = startpos+angle_inc*(0:iterations-1); % angle = p1:iterations:p2; angle_range = angle_inc*(iterations-1);
detect_E = zeros(iterations,1);
detect_E_std = zeros(iterations,1);
x_slice = zeros(iterations,length(snap_x));
y_slice = zeros(iterations,length(snap_y));
x_pos = 0;
y_pos = 0;

% OLD
% % matfile experiments: way too slow... only if no other way!
% snaps_obj = matfile('snaps_obj.mat','Writable',true);
% snaps_obj.snaps = zeros(2832,4244,iterations, 'single');

% Preallocate for speed and ensure data types
snaps = zeros(2832,4244,iterations, 'uint16'); % single is workable, double is too much, uint16 is best
snaps_full = zeros(2832,4244,iterations,avg, 'uint16'); % single is workable, double is too much, uint16 is best
snap_i = zeros(2832,4244,avg, 'uint16');
sum_i = zeros(avg,1,'single');

%% Phase 1: read, set and wait
f_setPos(s,startpos) % move to starting pos

pause(wait) % wait until the motor gets there
p1 = f_readPos(s); % Store position (check)

t1_dt = datetime; % store time
disp('start '); disp(t1_dt)
%% Phase 2: iterate over angle range
for iterate = 1:iterations
    for i = 1:avg
        snap_i(:,:,i) = (getsnapshot(vid)); % double % convert to double and store snaps: WHY? NO NEED TO MAKE DOUBLE! only increases size.. 
        % It is a 3D matrix that contains a "avg" number of images
        %% Filtering below LEVEL:
        snap_i(snap_i < level) = 0; % bootstrapping command
        sum_i(i) = sum(sum(snap_i(snap_x,snap_y,i))); % Calculate energy in filtered snaps
        % 
    end
    % Average of the "avg" number of images on the 3rd dimension (so between the images)
    snap = mean(snap_i(:,:,:),3); %,'native'); % Let convert to double to
    % allow better precision, uint16 loses precision
    
    %% Save shots (turn off during real meas, too slow)
    if capt || save_capt
        figure('Position',[250 100 700 500] ); imagesc(snap); colorbar; title('CCD image')
        if save_capt
            saveas(gcf,['shots/shot_' num2str(angle(iterate)) '.png']) % Save image
        end
        figure('Position',[1000 100 700 500] ); imagesc(log(snap)); colorbar; title('LOG CCD image')
        if save_capt
            saveas(gcf,['shots/shot_' num2str(angle(iterate)) '_log.png']) % Save image
        end
        %pause(1);
    end


    %% Save snaps in 3D matrix => This caused everything to go to shit... because it was DOUBLE!
    % RAM overflowed + Virtual Memory became too large (>12 GB!) % SOLVED!!!
    snaps(:,:,iterate) = snap;
%     snaps_full(:,:,iterate,1:avg) = snap_i;
    
%     % Check storing to HDD immediately, to prevent memory problems
%     % ONLY IF NECESSARY, much slower!
% %     save('snaps.mat','angle','snaps','-v7.3'): standard way, does not allow 'appending'. matfile does!
%     snaps_obj.snaps(1:2832,1:4244,iterate) = snap; % this seems to work now...

    %% Average and STD
    detect_E(iterate) = mean(sum_i); % Mean
    detect_E_std(iterate) = std(sum_i); % Standard deviation

    %% move stage with chosen angle
    f_movePos(s,angle_inc);
    pause(wait_inc) % Wait between the motor steps
    close all
end
delete('counter.mat'); set(0,'DefaultFigureVisible','on'); 
p2 = f_readPos(s); % End position

%% Time calc
% MATLAB built in:
t2_dt = datetime;
disp('stop '); disp(t2_dt)
time = t2_dt - t1_dt;
disp('It took '); % datestr(time,'SS') ' seconds'])
disp(time)

%% Phase 3: processing
% Don't remove noise, already removed by filtering (bootstrap-type)
% detect_E = detect_E - noise; % Remove background noise (dark current)

% Create plot and save
b = 1;      e = 0;

% Compensate for stronger TE: FIX WITH POLARIZER!
if strcmp(polarization, 'TE'); power = f_laser_intensity(source_power)*1; else; power = f_laser_intensity(source_power); end
if strcmp(pinhole, '40um'); power = power/28; elseif strcmp(pinhole, '100um'); power = power/100; end
detect_E_plot = detect_E/power; % fixed for real maximum
detect_E_std_plot = detect_E_std/power; %/mean(detect_E);
figure; errorbar(angle(b:end-e), detect_E_plot(b:end-e), detect_E_std_plot(b:end-e), 'o-')
title('Summed intensity over camera')
xlabel('Angle (^o)')
ylabel('Normalized instensity')
axis([angle(1) angle(end) 0 1])
saveas(gcf,[add '.png'])

STD = num2str(mean(detect_E_std_plot)*100);
disp(['The standard deviation on the measurement is:  ' STD ' %'])
fprintf('\n')
 
%% Saving
% Save values
save('detect_E.mat','angle','detect_E','detect_E_plot','detect_E_std',...
     'detect_E_std_plot','noise','source_power','snap_x','snap_y')
save(['STD = ' STD '.txt'],'STD')

% Save snaps: HEAVY SHIT!!! not anymore :)
tic
save('snaps.mat','angle','snaps','-v7.3')
toc
% tic
% save('snaps_full.mat','angle','snaps_full','-v7.3') % Used to save all
%                                                snaps, but not anymore :)
% toc
% clearvars snaps_full

cd .. % Return to previous folder