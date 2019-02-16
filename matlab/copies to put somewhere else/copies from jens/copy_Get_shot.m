% Model: Ximea xiD MD120MU-SY
close all; 
clearvars -except s vid vid_device detect_E detect_E_std angle; %clc
clearvars snap*
Parameters; % Always executes internally where it is required
set(0,'DefaultFigureVisible','on'); % off or on

% avg is the number of frames averaged over
vid.FramesPerTrigger = 1; % Default frames to capture per trigger
% Parameters
% avg_shot
% bootstrap

tic
snap_x = ':'; %1300:1700;
snap_y = ':'; % 200:3500;

%% Take snapshot(s)
snap_i = zeros(2832,4244); % Camera resolution, variable initialization
sum_i = ones(1); % variable initialization
for i = 1:avg_shot
    snap_i(:,:,i) = double(getsnapshot(vid));
    if bootstrap
    snap_i(snap_i < bootstrap) = 0; % bootstrapping command (new in V7)
    end
    sum_i(i) = sum(sum(snap_i(snap_x,snap_y,i)));
end
snap = mean(snap_i(snap_x,snap_y,:),3);

% figure('Position',[400 100 700 500] ); imagesc(snap_i(:,:,1)); colorbar;
figure('Position',[250 100 700 500] ); imagesc(snap); colorbar; title('CCD image')
figure('Position',[1000 100 700 500] ); imagesc(log(snap)); colorbar; title('LOG CCD image')
saveas(gcf,'shot.png') % Save image

try 
    source_power; 
catch
    source_power = f_LD_c2p(source_current); % Calculates with the source
                                             % current
end
toc

if show_parameters == 1
    snap_mean = mean(sum_i);
    disp(['Snap mean: ' num2str(snap_mean)]);

    snap_std = std(sum_i);
    disp(['Snap std: ' num2str(snap_std)]);

    Intensity_min_noise = snap_mean - f_power_normalizer(snap_x,snap_y);
    disp(['Intensity minimum noise: ' num2str(Intensity_min_noise)]);

    Predicted_laser_intensity = f_laser_intensity(source_power);
    disp(['Predicted laser intensity: ' num2str(Predicted_laser_intensity)]);
end

