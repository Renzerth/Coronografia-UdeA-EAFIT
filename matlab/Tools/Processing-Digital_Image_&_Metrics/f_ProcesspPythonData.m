function [] = f_ProcesspPythonData()

% Inputs:
%
%
% Outputs:

% Post-processing of the data (application of the metric of the degree of
% extintion)


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% INITIALIZATION
%% Processing initialization
infoDelim = '/';
fileFormat = '.mat';
[foundFiles] = loadFilesFromDir(fileFormat);
%%% Loading all the measurements
% 4 variables are loaded in the structure:

%%% Experiment information
tcvect = foundFiles.TCRanges;
glvect = foundFiles.GLRanges;

totalGL = length(glvect);
totalTC = length(tcvect);
totalImgs = totalGL*totalTC;
measInfo = cell(1,totalImgs);

%%% Processing Handlers
getIntensity = @(complexField) abs(fftshift(complexField)).^2;
scaleIntensity = @(intensity,maxVal) f_ScaleMatrixData(intensity,0,max(intensity(:))/maxVal);

%%%  Loading the reference measurement
refMeas = getIntensity(foundFiles.PSFreference);
maxRefVal = max(refMeas(:));
refMeas = f_ScaleMatrixData(refMeas,0,1);

%%% 'Experimental' images
expMeas = arrayfun(@(dataIndex) scaleIntensity(getIntensity(foundFiles.PSFoutputFields(:,:,dataIndex)),maxRefVal), 1:totalImgs,'UniformOutput',false);

%%% Data tags
idxgral = 1; % Init of general index that runs % on the range: [1,totalImgs]
for idxgl = 1:totalGL
    for idxtc = 1: totalTC
        tcstr = strcat('tc-',num2str(tcvect(idxtc)));
        glstr = strcat('ng-',num2str(glvect(idxgl)));
        measInfo{idxgral} = strcat(tcstr,'-',glstr); % Dataname for each
        idxgral = idxgral + 1; % The general index increases
    end
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Coordinates: Lambda/D, pixels and arcseconds
%% Find the center of the PSF image (with a binarization)
[~,~,aproxCenter,aproxRadius] = f_approximateSpotSize(refMeas);

%% Read the approximate center of the PSF reference
midX = aproxCenter(2);
midY = aproxCenter(1);

%% Cartesian coordinates with pixel units
% ORIGINAL
[ySize, xSize] = size(refMeas); % All images assumed of the same size as
% the refmeas
halfX = xSize/2;
halfY = ySize/2;

%% Centroid coordinates
% Centering shifting to the spot location
centerShiftX = (midX-halfX);
centerShiftY = (midY-halfY);

% Coordinates' origin set to the spot's center % 
xpixcenterd = (-halfX:halfX-1) - (centerShiftX - 1); % The x center is shifted 1
ypixcenterd = (-halfY:halfY-1) - (centerShiftY - 1); % The y center is shifted 1

%% Lambda over D scaling with the experimental spot size (THIS ONE'S USED)
% Pixel's size is scalled to the experimental spot's size
x = xpixcenterd/(aproxRadius);
y = ypixcenterd/(aproxRadius);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Profiles for the metrics
%% Reference Profile
%%% Metric-specific default parameters for the profile
oneSideProfile = 1; % Specifically needed for all the metrics. Ref: 1
shiftCart = [0,0]; % midX,midY already account for the shift

%%% Find the reference profile
disp('Calculating Profiles...');
[xrefProx,yrefProf,HprofRef,VprofRef,~,~] = f_makeImageProfile(x,y,midX,...
    midY,refMeas, shiftCart, oneSideProfile);
flippedAproxCenter = fliplr(aproxCenter); % [X,Y] Format
[averRefProfile] = f_getAverageRadialProfile(refMeas,[ySize, xSize], ...
    flippedAproxCenter);

%% Profile of the measurements
Hprofmeas = cell(1,totalImgs);
Vprofmeas = cell(1,totalImgs);
averMeasProfile = cell(1,totalImgs);

for idxgral = 1:totalImgs
    [~,~,Hprofmeas{idxgral},Vprofmeas{idxgral},~,~] = ...
        f_makeImageProfile(x,y,midX,midY,expMeas{idxgral},shiftCart, ...
        oneSideProfile);
    [averMeasProfile{idxgral}] = f_getAverageRadialProfile(...
        expMeas{idxgral},[ySize, xSize],flippedAproxCenter);
end

%% Measurement profile choice
metricProfile = 3;
switch metricProfile
    case 1 % Vertical profile
        radialIntensityRef = VprofRef;
        radialIntensityMeas = Vprofmeas; % One-sided
        cartcoord = yrefProf;
        titprof = '(vertical profile)';
        
    case 2 % Horizontal profile
        radialIntensityRef = HprofRef;
        radialIntensityMeas = Hprofmeas; % One-sided
        cartcoord = xrefProx;
        titprof = '(horizontal profile)';
        
    case 3 % Radial averaged profile (ImageAnalyst)
        radialIntensityRef = averRefProfile;
        radialIntensityMeas = averMeasProfile;
        trimRange = 1:min(numel(xrefProx),numel(yrefProf));
        radialAverageDist = mean([xrefProx(trimRange);
            yrefProf(trimRange)]);
        cartcoord = interp1(trimRange,radialAverageDist,1: ...
            length(averRefProfile)); % Spatial size interpolation range
        aproxRadius = find(cartcoord==1); % Airy range (corresponding
        % pixel) of interpolated data;
        titprof = '(radial average profile)';
        
    otherwise
        error('"metricProfile" must be either 1 or 2');
end
disp('Done.');

%% Data Cropping for view range
rangeFactor = 2.5; % Ref: 2 Number of Airy disks (Lambda/D times from center)
croppedMeasData = cell(1,totalImgs);
croppedCoorVect = x(abs(x)<=rangeFactor);
[cropRange] = f_computePSFCropRange(rangeFactor,2*aproxRadius,aproxCenter);
for idxgral = 1:totalImgs
    [croppedMeasData{idxgral}] = f_cropPSFrange(expMeas{idxgral},cropRange);
end

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Metric application
% Here, all the metrics are calculated but on a later stage only some are
% plotted

%% Processsing of profiles -- Encircled Energy Factor metric
disp('Calculating Encircled Energy Factor...');
measEEFcurves = cell(1,totalImgs);
measNormIntensity = cell(1,totalImgs);
dynamicProfileTitle = cell(1,totalImgs);
[refEEFcurve, refNormIntensity] = f_calculateEEF(radialIntensityRef);

for idxgral = 1:totalImgs
    [measEEFcurves{idxgral}, measNormIntensity{idxgral}] = ...
        f_calculateEEF(radialIntensityMeas{idxgral});
    dynamicProfileTitle{idxgral} = sprintf('of %s: %s',titprof, ...
        measInfo{idxgral});
end
disp('Done.');

%% Processsing of profiles -- Logarithmic SNR & Logarithmic RMS
disp('Calculating Logarithmic SNR and RMS...');
logSNR = cell(1,totalImgs);
logRMS = zeros(1,totalImgs);

for idxgral = 1:totalImgs
    [logSNR{idxgral}] = f_calculateLogSNR(radialIntensityMeas{idxgral}, ...
        radialIntensityRef);
    [logRMS(idxgral)] = f_calculateLogRMS(logSNR{idxgral});
end
disp('Done.');

%% Processsing of profiles -- Throughput, Power Suppression and logRMS
disp('Calculating Throughput...');
powerSupr = zeros(totalTC,totalGL);
arrangedEEF = cell(totalTC,totalGL);
arrangedProfiles = cell(totalTC,totalGL);
arrangedLogRMS = zeros(totalTC,totalGL);

for tcIndx = 1:totalTC
    for glIndx = 1:totalGL
        idxgral = tcIndx + (glIndx - 1)*totalTC; % Reversed width/index
        powerSupr(tcIndx,glIndx) = measEEFcurves{idxgral}(aproxRadius);
        arrangedEEF{tcIndx,glIndx} = measEEFcurves{idxgral};
        arrangedProfiles{tcIndx,glIndx}  = radialIntensityMeas{idxgral};
        arrangedLogRMS(tcIndx,glIndx) = logRMS(idxgral);
    end
end
arrangedProfiles = arrangedProfiles';
arrangedLogRMS = arrangedLogRMS';
disp('Done.');

%% Processsing of profiles -- Attenuation Ratios
disp('Calculating Attenuation Ratios...');
attenuationRatio = cell(1,totalImgs);
EEFattenuationRatio = cell(1,totalImgs);

for idxgral = 1:totalImgs
    [~,attenuationRatio{idxgral}] = f_calculateAttenuat( ...
        radialIntensityMeas{idxgral},radialIntensityRef);
    [~,EEFattenuationRatio{idxgral}] = f_calculateEEFAttenuat( ...
        measEEFcurves{idxgral},refEEFcurve);
end
disp('Done.');

%% Intensity profile gradient
disp('Calculating Intensity Gradient...');
measIntensityGradient = cell(1,totalImgs);
refIntensityGradient = gradient(radialIntensityRef) ;

for idxgral = 1:totalImgs
    measIntensityGradient{idxgral} = gradient( ...
        radialIntensityMeas{idxgral}); % gradient [returns n elements] or
    % diff [returns n-1 elements]
end
disp('Done.');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Metric Plotting
%% Plot Settings
xlab = 'Angular separation [\lambda/D]';
ylab = 'Angular separation [\lambda/D]';
tol = 0; % 0: no need to symmetrically truncate the profile. Ref: 0
plotData = 0; % Shows the profile lines. Ref: 1
plotH = 1;
plotV = 0;
metricSel = 3; % Type of metric -- BYPASS VARIABLE
                        % 1: Profiles
                        % 2: EEF: Encircled Energy Factor
                        % 3: Throughput (arranged EEF)
                        % 4: Power Supression
                        % 5: Relative contrast
                        % 6: Relative contrast (arranged)
                        % 7: Logarithmic SNR
                        % 8: Logarithmic RMS (RMS of logarithmic SNR) [arranged]
                        % 9: Attenuation Ratios
                        % 10: Gradient of Intensity
                        % 11: Plot Cropped intensity
                        % 12: Plot Images Mosaic

fontSize = 14; %[pts] Ref: 14
lineWidth = 1.5; %[pts] Ref: 1.5
colorSet = [1 0 0 ; 0 1 0; 0.8500 0.3250 0.0980; ...
    0 0 1; 0.9290 0.6940 0.1250; 0 1 1; 0.4940 0.1840 0.5560; ...
    1 0 1; 0.6350 0.0780 0.1840; 0 0 0];
lineStyle = '-.';
markerSet = [{'o'},{'+'},{'s'},{'>'},{'d'},{'x'},{'p'},{'^'},{'h'},{'v'}]';
plotSpec = arrayfun(@ (index) strcat(markerSet{index},lineStyle), ...
    1:length(markerSet),'UniformOutput',false); % Joints the line specs strings

%% Plot Selection
switch metricSel
    case 1
        %% Analysis Figures Plotting -- Profile Lines
        if metricProfile == 3 && plotData ~= 0
            warning(['Metric Profile 3 does not have a spatial related' ...
                'profile. Profile line disabled.']);
            plotData = 0;
        end
        disp('Plotting Profile Lines...');
        titRef = 'Profile of Reference Intensity';
        [refPoints] = getPlotCenterCoor([ySize, xSize] ,midX,midY,shiftCart);
        f_plotLinearProfiles(refMeas,x,y,cartcoord,cartcoord,titRef,xlab,ylab,plotData, radialIntensityRef, plotH,plotV,tol,refPoints,fontSize,lineWidth);
        
        for idxgral = 1:totalImgs
            f_plotLinearProfiles(expMeas{idxgral},x,y,cartcoord,cartcoord,dynamicProfileTitle{idxgral},xlab,ylab,plotData,radialIntensityMeas{idxgral},plotH,plotV,tol,refPoints,fontSize,lineWidth);
            fprintf('Plotting... %d/%d\n\r', idxgral, totalImgs);
        end
        
    case 2
        %% Analysis Figures Plotting -- Encircled Energy Factor metric
        disp('Plotting Encircled Energy Factor...');
        titprof = strcat(titprof,':',{' '},'Non-Coronagraphic');
        f_plotEEF(cartcoord,refEEFcurve,refNormIntensity,titprof,xlab, ...
            fontSize,lineWidth,colorSet,lineStyle,markerSet); % Reference
        
        for idxgral = 1:totalImgs
            f_plotEEF(cartcoord,measEEFcurves{idxgral}, ...
                measNormIntensity{idxgral},dynamicProfileTitle{idxgral},xlab, ...
                fontSize,lineWidth,colorSet,lineStyle,markerSet);
            fprintf('Plotting... %d/%d\n\r', idxgral, totalImgs);
        end
        
    case 3
        %% Analysis Figures Plotting -- Throughput
        disp('Plotting Throughput...');
        plotRange = 1:totalTC;
        plotAlotFunc = @(reference, Data, plotSpec, color,lineWidth) plot(reference, Data, plotSpec,'color', color, 'LineWidth', lineWidth);
        legendCell = cellstr(num2str(tcvect(plotRange)', 'TC=%d'));
        for indexGL = 1:totalGL
            figure('color', 'white');
            hold on; arrayfun(@(indexTC) plotAlotFunc(cartcoord,arrangedEEF{indexTC,indexGL},plotSpec{indexTC}, colorSet(indexTC,:),lineWidth),plotRange); hold off;
            xlabel('Angular separation (\lambda/D)','FontSize',fontSize,'FontWeight','bold');
            ylabel('Throughput (EEF)','FontSize',fontSize,'FontWeight','bold');  %  xlabel('Radial Distance (\lambda/D)')
            title(sprintf('Throughput of topological charges at NG = %d',glvect(indexGL)));
            set(gca,'FontSize',fontSize,'FontWeight','normal'); legend(legendCell,'Location','southeast'); grid on; axis square;
            fprintf('Plotting group... %d/%d\n\r', indexGL, totalGL); %xlim([0,2])
        end
        
    case 4
        %% Analysis Figures Plotting -- Power Supression in the Airy disk
        disp('Plotting Power Supression...');
        plotRange = 1:totalTC;
        figure('color', 'white');
        hold on; arrayfun(@(index) plot(glvect,powerSupr(index,:),plotSpec{index},'color',colorSet(index,:),'LineWidth',lineWidth), plotRange); hold off;
        axis fill; % They used to be too squared!
        title('Power supression of TCs at different phase levels','FontSize',fontSize,'FontWeight','bold');
        xlabel('Discretization level','FontSize',fontSize,'FontWeight','bold');
        ylabel('EEF in Airy disk','FontSize',fontSize,'FontWeight','bold'); % OLD:  ylabel('EEF at Airy Range');
        legendCell = cellstr(num2str(tcvect(plotRange)', 'TC=%d')); legend(legendCell); grid on; axis square;
        set(gca,'FontSize',fontSize,'FontWeight','normal')
        
    case 5
        %% Analysis Figures Plotting -- Relative Contrast
        disp('Plotting Relative Contrast...');
        for idxgral = 1:totalImgs
            f_plotContrast(cartcoord,radialIntensityRef,radialIntensityMeas{idxgral},dynamicProfileTitle{idxgral},fontSize,lineWidth)
            fprintf('Plotting... %d/%d\n\r', idxgral, totalImgs); %xlim([0,2.5]); % 2 Airy disks ylim([1e-3 1]); % Maximum attenuation
        end
        
    case 6
        %% Analysis Figures Plotting -- Relative Contrast (Arranged) [grouped gl's]
        disp('Plotting Arranged Relative Contrast...');
        plotRange = 1:totalGL;
        plotAlotFunc = @(reference, Data, plotSpec, color,lineWidth) plot(reference,Data,plotSpec,'color', color, 'LineWidth',lineWidth);
        legendCell = cellstr(num2str(glvect(plotRange)', 'Coronagraphic: NG=%d')); legendCell = [{'Non-Coronagraphic'}; legendCell];
        
        for indexTC = 1:totalTC
            figure('color', 'white');
            hold on; plot(cartcoord, radialIntensityRef,plotSpec{1},'color',colorSet(1,:),'LineWidth',lineWidth); set(gca,'yscale','log');
            arrayfun(@(indexGL) plotAlotFunc(cartcoord, arrangedProfiles{indexGL,indexTC},plotSpec{indexGL+1},colorSet(indexGL+1,:),lineWidth),plotRange); hold off;
            xlabel('Angular separation [\lambda/D]','FontSize',fontSize,'FontWeight','bold');
            ylabel('Relative contrast of the radial intensities [logscale]','FontSize',fontSize,'FontWeight','bold');
            title(sprintf('Raw Contrast NG Comparison with TC = %d',tcvect(indexTC)),'FontSize',fontSize,'FontWeight','bold');
            set(gca,'FontSize',fontSize,'FontWeight','normal'); legend(legendCell); grid on;
            fprintf('Plotting group... %d/%d\n\r', indexTC, totalTC); set(gca,'yscale','log');
%             xlim([0,2.5]); % 2 Airy disks
%             ylim([1e-3 1]); % Maximum attenuation
        end
        
    case 7
        %% Analysis Figures Plotting -- Logarithmic Signal-to-Noise Ratio
        disp('Plotting Logarithmic SNR...');
        for idxgral = 1:totalImgs
            f_plotLogSNR(cartcoord,logSNR{idxgral},dynamicProfileTitle{idxgral},xlab,fontSize,lineWidth)
            fprintf('Plotting... %d/%d\n\r', idxgral, totalImgs);
        end
        
    case 8
        %% Analysis Figures Plotting -- Logarithmic RMS
        disp('Plotting Logarithmic RMS...');
        plotRange = 1:totalGL;
        legendCell = cellstr(num2str(glvect(plotRange)', 'Coronagraphic: NG=%d'));
        
        figure('color', 'white');
        hold on; arrayfun(@(indexGL) plot(tcvect, arrangedLogRMS(indexGL,:), plotSpec{indexGL},'color',colorSet(indexGL,:),'LineWidth',lineWidth),plotRange);hold off
        title('Coronographic RMS analysis for GL effects','FontSize',fontSize,'FontWeight','bold');
        xlabel('Vortex Topological Charge','FontSize',fontSize,'FontWeight','bold');
        ylabel('Root Mean Square of Logarithmic SNR','FontSize',fontSize,'FontWeight','bold');
        legend(legendCell); set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on;
        
    case 9
        %% Analysis Figures Plotting -- Attenuation Ratios
        disp('Plotting Attenuation Ratios...');
        for idxgral = 1:totalImgs
            f_plotAttenRatios(cartcoord,attenuationRatio{idxgral},EEFattenuationRatio{idxgral},measInfo{idxgral},fontSize,lineWidth)
            fprintf('Plotting... %d/%d\n\r', idxgral, totalImgs);
        end
        
    case 10
        %%  Analysis Figures Plotting -- Gradient of Intensity
        disp('Plotting Gradient of Intensity...');
        for idxgral = 1:totalImgs
            f_plotGradient(cartcoord,measIntensityGradient{idxgral},measNormIntensity{idxgral},fontSize,lineWidth)
            fprintf('Plotting group... %d/%d\n\r', idxgral, totalImgs);
        end
        
    case 11
        %%  Analysis Figures Plotting -- Plot Cropped intensity
        disp('Plotting Cropped Images...');
        for idxgral = 1:totalImgs
            figure('color','white');
            imagesc(croppedCoorVect,croppedCoorVect,croppedMeasData{idxgral});
            xlabel('Angular separation (\lambda/D)','FontSize',fontSize,'FontWeight','bold');
            ylabel('Angular separation (\lambda/D)','FontSize',fontSize,'FontWeight','bold');
            title(sprintf('Coronagraphic PSF: %s',measInfo{idxgral}));
            set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on; axis square;
            colormap(viridis); colorbar; set(gca,'GridColor',[1,1,1]);
            fprintf('Plotting... %d/%d\n\r', idxgral, totalImgs);
        end
        
    case 12
        %%  Analysis Figures Plotting -- Plot Images Mosaic
        disp('Plotting Images Mosaic...');
        saveEnabled = false;
        titleSet = arrayfun(@(index) sprintf('TC:%d',tcvect(index)),1:totalTC,'UniformOutput',false);
        yLabelSet = arrayfun(@(index) sprintf('GL:%d',glvect(index)),1:totalGL,'UniformOutput',false);
        xLabelSet = cell(tcIndx,glIndx);
        
        arrangedCroppedImages = cell(totalTC,totalGL);
        for tcIndx = 1:totalTC
            for glIndx = 1:totalGL
                idxgral = tcIndx + (glIndx - 1)*totalTC; % Reversed width/index
                arrangedCroppedImages{tcIndx,glIndx}  = croppedMeasData{idxgral};
            end
        end
        
        f_plotMosaic(arrangedCroppedImages,croppedCoorVect,croppedCoorVect,titleSet,xLabelSet,yLabelSet,viridis,fontSize,saveEnabled)
        
    case 13 % WON'T BE USED
        %% Analysis Figures Plotting -- Mean Squared Error
        disp('Plotting Mean Squared Error...');
        fprintf('Underconstruction... %d/%d\n\r', 0, 0);
        % tit = 'Mean Squared Error';
        % % 13: MSE: Mean Squared Error
        
    otherwise
        warning('Unavailable Plot.');
end
disp('Rendering...');


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%  Save and finish the processing
%  This is not correct YET

%% Saving
%
% processedImgfullpath = strcat(processedImgname,measInfo{idxgral});
% % Explanation: saveas(variable,directory+filename,extension)
% saveas(gcf,strcat(processedImgfullpath),imgformat); % Saves the last shown figure

% OLD:
%   % Explanation: imwrite(variables,directory+filename+extension)
%   imwrite(expMeas{idxgral}, strcat(processedImgfullpath,dataformat));

%% End notification
N = 4; % Number of beeps for the processing
beepSound = 1;
f_EndBeeps(N,beepSound);

%% Ask to leave figures open
answer = questdlg('Do you want to close all the figures?','Processing finished','yes','no','no');
if strcmp(answer,'yes') % Compare string
    close all;
end

end