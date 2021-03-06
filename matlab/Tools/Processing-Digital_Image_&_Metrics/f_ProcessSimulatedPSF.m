function [] = f_ProcessSimulatedPSF()
% Inputs:

% Outputs:

% Post-processing of the data (application of the metric of the degree of extintion)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% INITIALIZATION
%% Processing initialization
fileFormat = '.mat';
[foundFiles] = f_loadFilesFromDir(fileFormat);
%%% Loading all the measurements
% 4 variables are loaded in the structure:

%%% Experiment information
tcvect = 1:foundFiles{1}.TC;
glvect = foundFiles{1}.glvect;

totalGL = length(glvect);
totalTC = length(tcvect);
totalImgs = totalGL*totalTC;
measInfo = cell(1,totalImgs);

%%% Processing Handlers
scaleIntensity = @(intensity,maxVal) f_ScaleMatrixData(intensity,0,max(intensity(:))/maxVal);

%%%  Loading the reference measurement
refMeas = foundFiles{1}.PSFreference;
maxRefVal = max(refMeas(:));
refMeas = f_ScaleMatrixData(refMeas,0,1);

%%% Simulated intensities
expMeas = cellfun(@(intensities) scaleIntensity(intensities,maxRefVal), foundFiles{1}.PSFplaneIntensities,'UniformOutput', false);

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
[ySize, xSize] = size(refMeas); % All images assumed of the same size as
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
[xrefProx,yrefProf,~,~,~,~] = f_makeImageProfile(x,y,midX,midY,refMeas, shiftCart, oneSideProfile);
flippedAproxCenter = fliplr(aproxCenter); % [X,Y] Format
[radialIntensityRef] = f_getAverageRadialProfile(refMeas,[ySize, xSize],flippedAproxCenter);

%% Profile of the measurements
radialIntensityMeas = cell(1,totalImgs);

for idxgral = 1:totalImgs
    [radialIntensityMeas{idxgral}] = f_getAverageRadialProfile(...
        expMeas{idxgral},[ySize, xSize],flippedAproxCenter);
    radialIntensityMeas{idxgral}(isnan(radialIntensityMeas{idxgral})) = 0;
end

%% Measurement profile choice Radial averaged profile 
trimRange = 1:min(numel(xrefProx),numel(yrefProf));
radialAverageDist = mean([xrefProx(trimRange); yrefProf(trimRange)]);
cartcoord = interp1(trimRange,radialAverageDist,1:length(radialIntensityRef)); % Spatial size interpolation range
aproxRadius = find(cartcoord==1); % Airy range (corresponding pixel) of interpolated data;
titleProf = '(radial average profile)';
disp('Done.');

%% Data Cropping for view range
rangeFactor = 1.5; % Ref: 2 Number of Airy disks (Lambda/D times from center)
croppedMeasData = cell(1,totalImgs);
croppedCoorVect = x(abs(x)<=2*rangeFactor); % Symmetric dual span coordinates cropping
[cropRange] = f_computePSFCropRange(rangeFactor,2*aproxRadius,aproxCenter);
for idxgral = 1:totalImgs
    [croppedMeasData{idxgral}] = f_cropPSFrange(expMeas{idxgral},cropRange);
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Metric application
% Here, all the metrics are calculated but on a later stage only some are plotted

%% Processsing of profiles -- Encircled Energy Factor metric
disp('Calculating Encircled Energy Factor...');
measEEFcurves = cell(1,totalImgs);
measNormIntensity = cell(1,totalImgs);
dynamicProfileTitle = cell(1,totalImgs);
[refEEFcurve, refNormIntensity] = f_calculateEEF(radialIntensityRef);

for idxgral = 1:totalImgs
    [measEEFcurves{idxgral}, measNormIntensity{idxgral}] = ...
        f_calculateEEF(radialIntensityMeas{idxgral});
    dynamicProfileTitle{idxgral} = sprintf('of %s: %s',titleProf, ...
        measInfo{idxgral});
end
disp('Done.');

%% Processsing of profiles -- Logarithmic SNR & Logarithmic RMS
disp('Calculating Logarithmic SNR and RMS...');
logSNR = cell(1,totalImgs);
logRMS = zeros(1,totalImgs);

for idxgral = 1:totalImgs
    [logSNR{idxgral}] = f_calculateLogSNR(radialIntensityMeas{idxgral}(1:aproxRadius), ...
        radialIntensityRef(1:aproxRadius));
    [logRMS(idxgral)] = f_calculateLogRMS(logSNR{idxgral});
end
disp('Done.');

%% Processsing of profiles -- Throughput, Power Suppression and logRMS
disp('Arranging Data Structures...');
powerSupr = zeros(totalTC,totalGL);
arrangedEEF = cell(totalTC,totalGL);
arrangedProfiles = cell(totalTC,totalGL);
arrangedLogSNR = cell(totalTC,totalGL);
arrangedLogRMS = zeros(totalTC,totalGL);

for tcIndx = 1:totalTC
    for glIndx = 1:totalGL
        idxgral = glIndx + (tcIndx - 1)*totalGL; % Reversed width/index
        powerSupr(tcIndx,glIndx) = measEEFcurves{idxgral}(aproxRadius);
        arrangedEEF{tcIndx,glIndx} = measEEFcurves{idxgral};
        arrangedProfiles{tcIndx,glIndx}  = radialIntensityMeas{idxgral};
        arrangedLogSNR{tcIndx,glIndx} = logSNR{idxgral};
        arrangedLogRMS(tcIndx,glIndx) = logRMS(idxgral);
    end
end
arrangedProfiles = arrangedProfiles';
arrangedLogSNR = arrangedLogSNR';
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
% refIntensityGradient = gradient(radialIntensityRef) ;

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
enableAxis = true;
metricSel = 13;

fontSize = 24; %[pts] Ref: 14
lineWidth = 1.5; %[pts] Ref: 1.5
markerSize = 5;
colorSet = [1 0 0 ; 0 1 0; 0.8500 0.3250 0.0980; ...
    0 0 1; 0.9290 0.6940 0.1250; 0 1 1; 0.4940 0.1840 0.5560; ...
    1 0 1; 0.6350 0.0780 0.1840; 0 0 0; 0.7216 0.1725 0.8784]; % lightBlue 0.3686 0.6941 0.7961 || dark green 0.1373 0.3255 0.2784
lineStyle = '-.';
xLimRange = [0,3];
yLimRange = [0,1];
markerSet = [{'o'},{'+'},{'s'},{'>'},{'d'},{'x'},{'p'},{'^'},{'h'},{'v'},{'<'}]';
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
        titleProf = strcat(titleProf,':',{' '},'Non-Coronagraphic');
        f_plotEEF(cartcoord,refEEFcurve,refNormIntensity,titleProf,xlab, ...
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
        plotAlotFunc = @(reference, Data, plotSpec, color,lineWidth) plot(reference, Data, plotSpec,'color', color, 'LineWidth', lineWidth,'MarkerSize',markerSize);
        legendCell = strtrim(cellstr(num2str(tcvect(plotRange)', 'TC=%d')));
        for indexGL = 1:totalGL
            figure('color', 'white');
            hold on; arrayfun(@(indexTC) plotAlotFunc(cartcoord,arrangedEEF{indexTC,indexGL},plotSpec{indexTC}, colorSet(indexTC,:),lineWidth),plotRange); hold off;
            xlabel('Angular separation (\lambda/D)','FontSize',fontSize,'FontWeight','bold');
            ylabel('Throughput (EEF)','FontSize',fontSize,'FontWeight','bold');  %  xlabel('Radial Distance (\lambda/D)')
            title(sprintf('Throughput of topological charges at GL = %d',glvect(indexGL)));
            set(gca,'FontSize',fontSize,'FontWeight','normal'); legend(legendCell,'Location','southeast'); grid on; axis square;
            fprintf('Plotting group... %d/%d\n\r', indexGL, totalGL); xlim(xLimRange);
            saveFigure(gcf, gca, [0,0,20,20], sprintf('Lyot_TC_throughput_GL_%d',glvect(indexGL)), 'svg');
        end
        
    case 4
        %% Analysis Figures Plotting -- Power Supression in the Airy disk
        disp('Plotting Power Supression...');
        plotRange = 1:totalTC;
        figure('color', 'white');
        hold on; arrayfun(@(index) plot(glvect,powerSupr(index,:),plotSpec{index},'color',colorSet(index,:),'LineWidth',lineWidth), plotRange); hold off;
        % axis fill; % They used to be too squared!
        title('Power supression of TCs at different phase levels','FontSize',fontSize,'FontWeight','bold');
        xlabel('Discretization level (GL)','FontSize',fontSize,'FontWeight','bold');
        ylabel('EEF in Airy disk','FontSize',fontSize,'FontWeight','bold'); % OLD:  ylabel('EEF at Airy Range');
        legendCell = strtrim(cellstr(num2str(tcvect(plotRange)', 'TC=%d')));
%         legend(legendCell); 
        columnlegend(2,legendCell,'FontSize',20,'Location','southwest');
        set(gca,'FontSize',fontSize,'FontWeight','normal'); xlim([12,256]); grid on; axis square;
%         saveFigure(gcf, gca, [0,0,20,20], 'EEF_power_suppression_Gaussian_N_3', 'svg');
        
    case 5
        %% Analysis Figures Plotting -- Relative Contrast
        disp('Plotting Relative Contrast...');
        for idxgral = 1:totalImgs
            f_plotContrast(cartcoord,radialIntensityRef,radialIntensityMeas{idxgral},dynamicProfileTitle{idxgral},fontSize,lineWidth)
            xlim(xLimRange); ylim(yLimRange);
            fprintf('Plotting... %d/%d\n\r', idxgral, totalImgs);
        end
        
    case 6
        %% Analysis Figures Plotting -- Relative Contrast (Arranged) [grouped gl's]
        disp('Plotting Arranged Relative Contrast...');
        plotRange = 1:totalGL;
        plotAlotFunc = @(reference, Data, plotSpec, color,lineWidth) plot(reference,Data,plotSpec,'color', color, 'LineWidth',lineWidth,'MarkerSize',markerSize);
        legendCell = strtrim(cellstr(num2str(glvect(plotRange)', 'GL=%d'))); legendCell = [{'GL=0'}; legendCell];
        
        for indexTC = 1:totalTC
            figure('color', 'white');
            hold on; plot(cartcoord, radialIntensityRef,plotSpec{1},'color',colorSet(1,:),'LineWidth',lineWidth,'MarkerSize',markerSize); set(gca,'yscale','log');
            arrayfun(@(indexGL) plotAlotFunc(cartcoord, arrangedProfiles{indexGL,indexTC},plotSpec{indexGL+1},colorSet(indexGL+1,:),lineWidth),plotRange); hold off;
            xlabel('Angular separation [\lambda/D]','FontSize',fontSize,'FontWeight','bold');
            ylabel('Logarithmic Intensity Contrast','FontSize',fontSize,'FontWeight','bold');
            title(sprintf('Radial PSF with Discretized Masks with TC = %d',tcvect(indexTC)),'FontSize',fontSize,'FontWeight','bold');
            set(gca,'FontSize',fontSize,'FontWeight','normal'); legend(legendCell); grid on; axis square;
            fprintf('Plotting group... %d/%d\n\r', indexTC, totalTC); set(gca,'yscale','log');
            xlim(xLimRange); % 2 Airy disks
            ylim(yLimRange); % Maximum attenuation
            saveFigure(gcf, gca, [0,0,20,20], sprintf('Intensity_Contrast_GL_high_TC_%d',indexTC), 'svg');
        end
        
    case 7
        %% Analysis Figures Plotting -- Logarithmic Signal-to-Noise Ratio
        disp('Plotting Logarithmic SNR...');
        for idxgral = 1:totalImgs
            f_plotLogSNR(cartcoord,logSNR{idxgral},dynamicProfileTitle{idxgral},xlab,fontSize,lineWidth)
            fprintf('Plotting... %d/%d\n\r', idxgral, totalImgs);
        end
        
    case 8
        %% Analysis Figures Plotting -- logSNR (Arranged) [grouped gl's]
        disp('Plotting Logarithmic SNR...');
        plotRange = 1:totalGL;
        plotAlotFunc = @(reference, Data, plotSpec, color,lineWidth) plot(reference,Data,plotSpec,'color', color, 'LineWidth',lineWidth,'MarkerSize',markerSize);
        legendCell = strtrim(cellstr(num2str(glvect(plotRange)', 'GL=%d')));
        
        for indexTC = 1:totalTC
            figure('color', 'white');
            hold on; arrayfun(@(indexGL) plotAlotFunc(cartcoord, arrangedLogSNR{indexGL,indexTC},plotSpec{indexGL+1},colorSet(indexGL+1,:),lineWidth),plotRange); hold off;
            xlabel('Log Scale Angular separation [\lambda/D]','FontSize',fontSize,'FontWeight','bold');
            ylabel('Logarithmic SNR','FontSize',fontSize,'FontWeight','bold');
            title(sprintf('Gray Level LSNR Comparison for TC = %d',tcvect(indexTC)),'FontSize',fontSize,'FontWeight','bold');
            set(gca,'FontSize',fontSize,'FontWeight','normal'); legend(legendCell,'Location','northeast'); grid on; axis square;
            fprintf('Plotting group... %d/%d\n\r', indexTC, totalTC); set(gca,'xscale','log');
            xlim([0,10]); %ylim([-4.5,-1]);
            saveFigure(gcf, gca, [0,0,20,20], sprintf('logSNR_TC_%d',indexTC), 'svg');
        end
        
    case 9
        %% Analysis Figures Plotting -- Logarithmic RMS
        disp('Plotting Logarithmic RMS...');
        plotRange = 1:totalGL;
        legendCell = strtrim(cellstr(num2str(glvect(plotRange)', 'GL=%d')));
        
        figure('color', 'white');
        hold on; arrayfun(@(indexGL) plot(tcvect, arrangedLogRMS(indexGL,:), plotSpec{indexGL},'color',colorSet(indexGL,:),'LineWidth',lineWidth),plotRange);hold off
        title('Coronagraphic RMS analysis for GL effects','FontSize',fontSize,'FontWeight','bold');
        xlabel('Vortex Topological Charge (TC)','FontSize',fontSize,'FontWeight','bold');
        ylabel('Root Mean Square of the Logarithmic SNR','FontSize',fontSize,'FontWeight','bold');
%         legend(legendCell,'Location','southeast'); 
        columnlegend(3,legendCell,'FontSize',24);
        set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on; axis square;
        xlim([tcvect(1),tcvect(end)]);
%         saveFigure(gcf, gca, [0,0,20,20], sprintf('Logarithmic_RMS_Gaussian_N_4%d',1), 'svg');
        
    case 10
        %% Analysis Figures Plotting -- Attenuation Ratios
        disp('Plotting Attenuation Ratios...');
        for idxgral = 1:totalImgs
            f_plotAttenRatios(cartcoord,attenuationRatio{idxgral},EEFattenuationRatio{idxgral},measInfo{idxgral},fontSize,lineWidth)
            fprintf('Plotting... %d/%d\n\r', idxgral, totalImgs);
        end
        
    case 11
        %%  Analysis Figures Plotting -- Gradient of Intensity
        disp('Plotting Gradient of Intensity...');
        for idxgral = 1:totalImgs
            f_plotGradient(cartcoord,measIntensityGradient{idxgral},measNormIntensity{idxgral},fontSize,lineWidth)
            fprintf('Plotting group... %d/%d\n\r', idxgral, totalImgs);
        end
        
    case 12
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
        
    case 13
        %%  Analysis Figures Plotting -- Plot Images Mosaic
        disp('Plotting Images Mosaic...');
        saveEnabled = false;
        logViewEnabled = true;
        fontSize = 17;
        titleSet = arrayfun(@(index) sprintf('TC:%d',tcvect(index)),1:totalTC,'UniformOutput',false);
        yLabelSet = arrayfun(@(index) sprintf('GL:%d',glvect(index)),1:totalGL,'UniformOutput',false);
        xLabelSet = cell(totalTC,totalGL);
        [croppedRefData] = f_cropPSFrange(refMeas,cropRange);
        
        arrangedCroppedImages = cell(totalTC,totalGL);
        for tcIndx = 1:totalTC
            for glIndx = 1:totalGL
                idxgral = glIndx + (tcIndx - 1)*totalGL; % Reversed width/index
                arrangedCroppedImages{tcIndx,glIndx}  = croppedMeasData{idxgral};
            end
        end
        
        if logViewEnabled == true
            logIntensity =  @(intensityData) 10*log10(intensityData);
            logCroppedRefData = logIntensity(croppedRefData);
            scalingLimits = computeScaleRange(logCroppedRefData,[0.2,1],[-15,0]);
            arrangedCroppedImages = cellfun(@(cellData) logIntensity(cellData), arrangedCroppedImages, 'UniformOutput', false);
        else
            scalingLimits = computeScaleRange(croppedRefData,[1,1]);
        end
        f_plotMosaic(arrangedCroppedImages,croppedCoorVect,croppedCoorVect,titleSet,xLabelSet,yLabelSet,viridis,fontSize,saveEnabled,enableAxis,scalingLimits)
        
    case 14
        %% Analysis Figures Plotting -- Gray level improvement 
        disp('Plotting Logarithmic RMS...');
        plotRange = 1:1:totalTC;
        legendCell = strtrim(cellstr(num2str(tcvect(plotRange)', 'TC=%d')));
        relativePercentage = @(matrixData) (matrixData(:,:) - repmat(matrixData(1,:), [size(matrixData,1),1]) )./repmat(matrixData(1,:), [size(matrixData,1),1])*100;
        GLlmprovement = relativePercentage(arrangedLogRMS);
        
        figure('color', 'white');
        hold on; arrayfun(@(indexTC) plot(glvect, GLlmprovement(:,indexTC), plotSpec{indexTC},'color',colorSet(indexTC,:),'LineWidth',lineWidth),plotRange);hold off
        title('Averaged Gray Level Improvement Effect','FontSize',fontSize,'FontWeight','bold');
        xlabel('Discretization level','FontSize',fontSize,'FontWeight','bold');
        ylabel('LRMS Improvement [%]','FontSize',fontSize,'FontWeight','bold');
%         legend(legendCell,'Location','northeast');
        columnlegend(2,legendCell,'FontSize',18,'Location','southeast');
        set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on; axis square; xlim([2,10])
%         saveFigure(gcf, gca, [0,0,20,20], sprintf('RMS_gray_level_improvement_%d',1), 'svg');
        
    case 15
        %% Cropped Referential PSF View -- LogView
        [croppedRefData] = f_cropPSFrange(refMeas,cropRange);
        
        figure('color','white');
        imagesc(croppedCoorVect,croppedCoorVect,croppedRefData);
        xlabel('Angular separation (\lambda/D)','FontSize',fontSize,'FontWeight','bold');
        ylabel('Angular separation (\lambda/D)','FontSize',fontSize,'FontWeight','bold');
        title('Referential Non-Coronagraphic PSF');
        set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on; axis square;
        colormap(viridis); colorbar; set(gca,'GridColor',[1,1,1]);
        saveFigure(gcf, gca, [0,0,20,20], 'Referential_PSF_normal_view', 'svg');
        
        figure('color','white');
        imagesc(croppedCoorVect,croppedCoorVect,10*log10(croppedRefData));
        xlabel('Angular separation (\lambda/D)','FontSize',fontSize,'FontWeight','bold');
        ylabel('Angular separation (\lambda/D)','FontSize',fontSize,'FontWeight','bold');
        title('Log View Referential Non-Coronagraphic PSF');
        set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on; axis square;
        colormap(viridis); colorbar; set(gca,'GridColor',[1,1,1]);
        saveFigure(gcf, gca, [0,0,20,20], 'Referential_PSF_log_view', 'svg');
        
    case 16 % WON'T BE USED
        %% Analysis Figures Plotting -- Mean Squared Error
        disp('Plotting Mean Squared Error...');
        fprintf('Underconstruction... %d/%d\n\r', 0, 0);
        % tit = 'Mean Squared Error';
        % % 13: MSE: Mean Squared Error
        
    otherwise
        warning('Unavailable Plot.');
end
disp('Rendering...');


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