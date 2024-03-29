function [x,y,radialIntensityRef,radialIntensityMeas,cartcoord, ...
refEEFcurve,refNormIntensity,measEEFcurves,measNormIntensity,logSNR, ...
attenuationRatio,EEFattenuationRatio,refIntensityGradient, ...
measIntensityGradient] = f_ProcessData(measfullpath,refmeasfullpath, ...
ProcessedDir,dataDir,pathSep,totalImgs,metricProfile, beepSound,L,NA...
,PP,measSimulated,glvect,tcvect)

% Inputs:
%
%
% Outputs:

% Post-processing of the data (application of the metric of the degree of
% extintion)



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% INITIALIZATION
%% Processing initialization

%%% Measure the processing time
% Copyright PhD student Jens de Pelsmaeker VUB B-PHOT 2018,Brussels,Belgium
t1_dt = datetime; % store time
disp('Processing started:'); disp(t1_dt)

% Saving on DEMAND.
savingEnabled = false;

% processedImgname = strcat(ProcessedDir,pathSep,'processed',infoDelim, ...
% cameraPlane,infoDelim); % Used to save the metrics                       (NOT USED RIGHT NOW)

%%% Loading all the measurements
% Explanation: load(directory+filename,variables)
struct = load(measfullpath); % Loads all the measured images & info
% Two variables are loaded in the structure:
% "expImgs" and "MeasInfo"

%%% Experimental images
expMeas = struct.expImgs;

%%% Experiment information
measInfo = struct.MeasInfo;

%%%  Loading the reference measurement
% The image is read in a uint8 format: integer with values that are
% normally in [0,255] (8-bit depth or dynamic range)
refMeas = imread(refmeasfullpath);
% since refMeas is a bmp image, it is loaded as uint8
    
%%% UINT8 format to Double for the reference image
% im2double duplicates the precision of the exponent leaving intact the
% mantisa. It as floating-point format that normalizes the images and 
% this operation is made on each RGB channel. rgb2gray does a similar 
% operation but scaling to a gray scale, where it  converts RGB values 
% to grayscale values by forming a weighted sum of the R, G, and B 
% components: 0.2989 * R + 0.5870 * G + 0.1140 * B
refMeas = im2double(refMeas);

% Lyotimg = refMeas;                                                           % TO BE USED FOR LYOT's METRICS
    
%% Debug example measurements (not used right now)
exampleDataDebug = 0;
if exampleDataDebug == 1
  %% Example images to process
  % Lyotimg = imread(strcat(dataDir,pathSep,'0_ExampleData',pathSep, ...
  %                         'data_ref_1.bmp')); % Lyot image
  % Lyotimg = rgb2gray(Lyotimg);
  refMeas = imread(strcat(dataDir,pathSep,'0_ExampleData',pathSep, ...
                                    'data_ref_2.bmp')); % PSF reference
  expMeas = {0,0}; % Cell initialization
  expMeas{1} = imread(strcat(dataDir,pathSep,'0_ExampleData', ...
                        pathSep,'data_ref_3.png')); % PSF measurement 1
  expMeas{2} = imread(strcat(dataDir,pathSep,'0_ExampleData', ...
                        pathSep,'data_ref_4.png')); % PSF measurement 2
  measInfo = {'data_ref_3','data_ref_3'};
  totalImgs = 2; % For the case of the simulated measurement
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Coordinates: Lambda/D, pixels and arcseconds
%% Find the center of the PSF image (with a binarization)
if measSimulated == 0
    Threshold = 0.3; % 70% of Airy's energy
    enableAxis = true;
else
    Threshold = 1; % Reference is a binary mask
    enableAxis = false;
end
[~,~,aproxCenter,aproxRadius] = f_approximateSpotSize(refMeas,Threshold);

%% Read the approximate center of the PSF reference
midX = aproxCenter(2); % Centroid of the reference in x
midY = aproxCenter(1); % Centroid of the reference in y
aproxCenter = [midY midX]; % shift of one saved back on aproxCenter

%% Cartesian coordinates with pixel units
% ORIGINAL
[ySize, xSize] = size(refMeas); % All images assumed of the same size as
% the refmeas
halfX = xSize/2; % Half size of the whole image in x
halfY = ySize/2; % Half size of the whole image in y

%% Centroid coordinates
% Centering shifting to the spot location
centerShiftX = (midX-halfX);
centerShiftY = (midY-halfY);

% Coordinates' origin set to the spot's center
xpixcenterd = (-halfX:halfX - 1) - (centerShiftX - 1); 
ypixcenterd = (-halfY:halfY - 1) - (centerShiftY - 1);
% The x,y centers are shifted 1 and this is substracted

%% Lambda over D scaling with the experimental spot size (THIS ONE'S USED)
% Pixel's size is scalled to the experimental spot's size
x  = xpixcenterd/(aproxRadius); 
y  = ypixcenterd/(aproxRadius);

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Profiles for the metrics
%% Reference Profile
%%% Metric-specific default parameters for the profile
oneSideProfile = 1; % Specifically needed for all the metrics. Ref: 1
shiftCart = [0,0]; % midX,midY already account for the shift. Ref: [0,0]

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
if measSimulated == 0
    rangeFactor = 0.5; % Ref: 1.5 Number of Airy disks (Lambda/D times from
                       % the center)
   Diameter = 2*aproxRadius;
else
    rangeFactor = 1; % The masks ocupy the whole screen
    Diameter = aproxRadius; % This is a "diameter"
end
croppedMeasData = cell(1,totalImgs); % same size as the measured bank
croppedCoorVect = x(abs(x)<=2*rangeFactor); % Symmetric dual span coordinates cropping

[cropRange] = f_computePSFCropRange(rangeFactor,Diameter,aproxCenter);
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
    measEEFcurves{idxgral}(isnan(measEEFcurves{idxgral})) = 0;
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

%% Processsing of profiles -- Throughput, Power Suppression, logSNR and logRMS
disp('Arranging Data Structures...');
totalGL = length(glvect);
totalTC = length(tcvect);
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
        arrangedProfiles{tcIndx,glIndx} = radialIntensityMeas{idxgral};
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
xlab = 'Radial Distance [a/\rho]';
ylab = 'Radial Distance [a/\rho]';
tol = 0; % 0: no need to symmetrically truncate the profile. Ref: 0
plotData = 0; % Shows the profile lines. Ref: 1
plotH = 1;
plotV = 0;
enableAxis = true;
metricSel = 13;
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
                
fontSize = 24; %[pts] Ref: 14
lineWidth = 1.3; %[pts] Ref: 1.5
markerSize = 4;
colorSet = [1 0 0 ; 0 1 0; 0.8500 0.3250 0.0980; ...
    0 0 1; 0.9290 0.6940 0.1250; 0 1 1; 0.4940 0.1840 0.5560; ...
    1 0 1; 0.6350 0.0780 0.1840; 0 0 0; 0.7216 0.1725 0.8784]; % lightBlue 0.3686 0.6941 0.7961 || dark green 0.1373 0.3255 0.2784
lineStyle = '-.';
xLimRange = [0,3];
yLimRange = [0,1.4e-3];
markerSet = [{'o'},{'+'},{'s'},{'>'},{'d'},{'x'},{'p'},{'^'},{'h'},{'v'},{'<'}]';
plotSpec = arrayfun(@ (index) strcat(markerSet{index},lineStyle), ...
    1:length(markerSet),'UniformOutput',false); % Joints the line specs strings

if measSimulated == 0
    colorM = viridis;
else
    colorM = gray;
end

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
        [refPoints] = f_getPlotCenterCoor([ySize, xSize] ,midX,midY, ...
                                           shiftCart);
        f_plotLinearProfiles(refMeas,x,y,cartcoord,cartcoord,titRef, ...
        xlab,ylab,plotData, radialIntensityRef, plotH,plotV,tol, ...
        refPoints,fontSize,lineWidth);
        
        for idxgral = 1:totalImgs
            f_plotLinearProfiles(expMeas{idxgral},x,y,cartcoord, ...
            cartcoord,dynamicProfileTitle{idxgral},xlab,ylab,plotData, ...
            radialIntensityMeas{idxgral},plotH,plotV,tol,refPoints, ...
            fontSize,lineWidth);
            fprintf('Plotting... %d/%d\n\r', idxgral, totalImgs);
        end
        
    case 2
        %% Analysis Figures Plotting -- Encircled Energy Factor metric
        disp('Plotting Encircled Energy Factor...');
        titprof = strcat(titprof,':',' Non-Coronagraphic');
        f_plotEEF(cartcoord,refEEFcurve,refNormIntensity,titprof,xlab, ...
            fontSize,lineWidth,colorSet,lineStyle,markerSet); % Reference
        
        for idxgral = 1:totalImgs
            f_plotEEF(cartcoord,measEEFcurves{idxgral}, ...
            measNormIntensity{idxgral},dynamicProfileTitle{idxgral}, ...
            xlab,fontSize,lineWidth,colorSet,lineStyle,markerSet);
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
            fprintf('Plotting group... %d/%d\n\r', indexGL, totalGL); xlim(xLimRange); % set(gca,'yscale','log');
            saveFigure(gcf, gca, [0,0,20,20], sprintf('PSF_LC2002_EXP_Throughput_GL_%d',glvect(indexGL)), 'svg');
        end
        
    case 4
        %% Analysis Figures Plotting -- Power Supression in the Airy disk
        disp('Plotting Power Supression...');
        plotRange = 1:totalTC;
        figure('color', 'white');
        powerSuprNoNan = powerSupr;
        powerSuprNoNan(isnan(powerSuprNoNan)) = 0;
        hold on; arrayfun(@(index) plot(glvect,powerSuprNoNan(index,:),plotSpec{index},'color',colorSet(index,:),'LineWidth',lineWidth,'MarkerSize',markerSize), plotRange); hold off;
        % axis fill; % They used to be too squared!
        title('Power supression of TCs at different phase levels','FontSize',fontSize,'FontWeight','bold');
        xlabel('Discretization level (GL)','FontSize',fontSize,'FontWeight','bold');
        ylabel('EEF in Airy disk','FontSize',fontSize,'FontWeight','bold'); % OLD:  ylabel('EEF at Airy Range');
        legendCell = strtrim(cellstr(num2str(tcvect(plotRange)', 'TC=%d'))); legend(legendCell,'Location','northeast'); grid on; axis square;
        set(gca,'FontSize',fontSize,'FontWeight','normal'); %xlim([2,12]);
        saveFigure(gcf, gca, [0,0,20,20], 'PSF_EXP_EEF_power_suppression_GL_high', 'svg');
        
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
            ylabel('Relative contrast of the radial intensities [logscale]','FontSize',fontSize,'FontWeight','bold');
            title(sprintf('Raw Contrast GL Comparison with TC = %d',tcvect(indexTC)),'FontSize',fontSize,'FontWeight','bold');
            set(gca,'FontSize',fontSize,'FontWeight','normal'); legend(legendCell); legend boxoff; grid on; axis square;
            fprintf('Plotting group... %d/%d\n\r', indexTC, totalTC); set(gca,'yscale','log');
            xlim(xLimRange); % 2 Airy disks % ylim(yLimRange); % Maximum attenuation
            saveFigure(gcf, gca, [0,0,20,20], sprintf('PSF_LC2002_EXP_intensities_GL_low_TC_%d',tcvect(indexTC)), 'svg');
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
            set(gca,'FontSize',fontSize,'FontWeight','normal'); legend(legendCell,'Location','southwest'); grid on; axis square;
            fprintf('Plotting group... %d/%d\n\r', indexTC, totalTC); set(gca,'xscale','log');
%             xlim(logSepRange); % 2 Airy disks
%             saveFigure(gcf, gca, [0,0,20,20], sprintf('figure_%d',indexTC), 'svg');
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
        legend(legendCell,'Location','northwest'); set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on; axis square;
        xlim([tcvect(1),tcvect(end)]);
        saveFigure(gcf, gca, [0,0,20,20], sprintf('PSF_EXP_Logarithmic_RMS_GL_high%d',1), 'svg');
        
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
            colormap(colorM); colorbar; set(gca,'GridColor',[1,1,1]);
            fprintf('Plotting... %d/%d\n\r', idxgral, totalImgs);
        end
        
    case 13
        %%  Analysis Figures Plotting -- Plot Images Mosaic
        disp('Plotting Images Mosaic...');
        saveEnabled = false;
        logViewEnabled = true;
        fontSize = 15;
        titleSet = arrayfun(@(index) sprintf('TC:%d',tcvect(index)),1:totalTC,'UniformOutput',false);
        yLabelSet = arrayfun(@(index) sprintf('GL:%d',glvect(index)),1:totalGL,'UniformOutput',false);
        xLabelSet = cell(totalTC,totalGL);
        [croppedRefData] = f_cropPSFrange(refMeas,cropRange);
        
        arrangedCroppedImages = cell(totalTC,totalGL);
        for tcIndx = 1:totalTC
            for glIndx = 1:totalGL
                idxgral = glIndx + (tcIndx - 1)*totalGL; % Reversed width/index
                arrangedCroppedImages{tcIndx,glIndx} = croppedMeasData{idxgral};
            end
        end
        
        if measSimulated == 0
            if logViewEnabled == true
                logIntensity =  @(intensityData) 10*log10(intensityData);
                logCroppedRefData = logIntensity(croppedRefData);
                scalingLimits = computeScaleRange(logCroppedRefData,[1,1],[0,0]);
                arrangedCroppedImages = cellfun(@(cellData) logIntensity(cellData), arrangedCroppedImages, 'UniformOutput', false);
            else
                scalingLimits = computeScaleRange(croppedRefData,[1,1]);
            end
            f_plotMosaic(arrangedCroppedImages,croppedCoorVect,croppedCoorVect,titleSet,xLabelSet,yLabelSet,colorM,fontSize,saveEnabled, enableAxis,scalingLimits) 
        else
            f_plotMosaic_angle(arrangedCroppedImages,croppedCoorVect,croppedCoorVect,titleSet,xLabelSet,yLabelSet,colorM,fontSize,saveEnabled,enableAxis) 
        end
        
    case 14
        %% Analysis Figures Plotting -- Logarithmic RMS GL Improvement
        disp('Plotting Logarithmic RMS...');
        plotRange = 1:1:totalTC;
        legendCell = strtrim(cellstr(num2str(tcvect(plotRange)', 'TC=%d')));
        %         relativePercentage = @(matrixData) (matrixData(:,:) - matrixData(1,:) )./matrixData(1,:)*100;
        relativePercentage = @(matrixData) (matrixData(:,:) - repmat(matrixData(1,:), [size(matrixData,1),1]) )./repmat(matrixData(1,:), [size(matrixData,1),1])*100;
        GLlmprovement = relativePercentage(arrangedLogRMS);

        figure('color', 'white');
        hold on; arrayfun(@(indexTC) plot(glvect, GLlmprovement(:,indexTC), plotSpec{indexTC},'color',colorSet(indexTC,:),'LineWidth',lineWidth),plotRange);hold off
        title('Averaged Gray Level Improvement Effect','FontSize',fontSize,'FontWeight','bold');
        xlabel('Discretization level','FontSize',fontSize,'FontWeight','bold');
        ylabel('LRMS Improvement [%]','FontSize',fontSize,'FontWeight','bold');
        legend(legendCell,'Location','northeast'); set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on;
%         saveFigure(gcf, gca, [0,0,20,20], sprintf('Logarithmic_RMS_%d',1), 'svg');

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
        saveFigure(gcf, gca, [0,0,20,20], 'Referential_LC2002_EXP_PSF_normal_view', 'svg');
        
        figure('color','white');
        imagesc(croppedCoorVect,croppedCoorVect,10*log10(croppedRefData));
        xlabel('Angular separation (\lambda/D)','FontSize',fontSize,'FontWeight','bold');
        ylabel('Angular separation (\lambda/D)','FontSize',fontSize,'FontWeight','bold');
        title('Log View Referential Non-Coronagraphic PSF');
        set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on; axis square;
        colormap(viridis); colorbar; set(gca,'GridColor',[1,1,1]);
        saveFigure(gcf, gca, [0,0,20,20], 'Referential_LC2002_EXP_PSF_log_view', 'svg');
    
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

if savingEnabled == true
    save(strcat(ProcessedDir,pathSep,'processedWorkspace.mat'));
end

%% End of the processing
% Author: PhD student Jens de Pelsmaeker VUB B-PHOT 2018, Brussels, Belgium
% MATLAB built in:
t2_dt = datetime;
disp('Processing finished:'); disp(t2_dt)
time = t2_dt - t1_dt;
disp('Processing took: '); % datestr(time,'SS') ' seconds'])
disp(time);

%% End notification
N = 4; % Number of beeps for the processing
f_EndBeeps(N,beepSound);

%% Ask to leave figures open
answer = questdlg('Do you want to close all the figures?','Processing finished','yes','no','no');
if strcmp(answer,'yes') % Compare string
    close all;
end


%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%% Airy size estimatives
%% Theoretical: diffraction-limited systems
% NA: numerical aperture of the lens
% L: wavelength in um
PP = PP*1e-6; % um to m % Measurement camera's PP

AiryDiskSpatial = 0.61*L/NA; % um
AiryDiskSpatialX = AiryDiskSpatial;
AiryDiskSpatialY = AiryDiskSpatial;
disp(strcat('Theoretical Airy"s radius in um:', {' '}, num2str(AiryDiskSpatialX), ' (x) and', {' '},  num2str(AiryDiskSpatialY), ' (y)'));

AiryDiskPixX = round(AiryDiskSpatialX/PP); % PP: camera's pixel pitch
AiryDiskPixY = round(AiryDiskSpatialY/PP);
disp(strcat('Theoretical Airy"s radius in pix:', {' '},  num2str(AiryDiskPixX), ' (x) and', {' '},  num2str(AiryDiskPixY), ' (y)'));

% In reality, there are aberrations and the real radius can be of
% about 60 times the theoretical one

%% Airy radius measured from the reference image
AiryDiskPixX = aproxRadius; % Just an example
AiryDiskPixY= aproxRadius; % Just an example
AiryDiskSpatialX = AiryDiskPixX*PP; % PP: camera's pixel pitch
AiryDiskSpatialY = AiryDiskPixY*PP;
disp(strcat('Estimated Airy"s radius in um:', {' '},  num2str(AiryDiskSpatialX), ' (x) and', {' '},  num2str(AiryDiskSpatialY), ' (y)'));
disp(strcat('Estimated Airy"s radius in pix:', {' '},  num2str(AiryDiskPixX), ' (x) and', {' '},  num2str(AiryDiskPixY), ' (y)'));

%% Airy radius from the EEF factor
% When it is the 70%
distances = cartcoord(refEEFcurve<=0.7);
AiryRadius = distances(end)*aproxRadius;
fprintf('Estimated EEF Airy"s radius in [Pixels]: %d.\n',AiryRadius);
% PENDIGL

end