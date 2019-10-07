%% Plot Settings
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
plotSpec = arrayfun(@ (index) strcat(markerSet{index},lineStyle),1:length(markerSet),'UniformOutput',false); % Joints the line specs strings

%% Plot ticks data
apertureDiameter = 3; % [mm]
gaussianWidth = apertureDiameter./[4, 3, 2, 1, 0.5, 0.2];

plotRangeA = 1:length(tcvect);
legendCellA = strtrim(cellstr(num2str(tcvect(plotRangeA)', 'TC=%d')));

plotRangeB = 1:length(gaussianWidth);
legendCellB = strtrim(cellstr(num2str(gaussianWidth(plotRangeB)', 'w=%1.2f')));

%% Plots -- Arranged LogRMS per Gaussian witdh per Gray Level
for indexGL = 1:length(glvect)
    perGLvariation = [AA(indexGL,:,1); AA(indexGL,:,2); AA(indexGL,:,3); AA(indexGL,:,4); AA(indexGL,:,5); AA(indexGL,:,6)]';
    figure('color', 'white');
    hold on; arrayfun(@(index) plot(gaussianWidth,perGLvariation(index,:),plotSpec{index},'color',colorSet(index,:),'LineWidth',lineWidth), plotRangeA); hold off;
    title(sprintf('Effect of Gaussian illumination-GL:%d',glvect(indexGL)),'FontSize',fontSize,'FontWeight','bold');
    xlabel('Gaussian Width [mm]','FontSize',fontSize,'FontWeight','bold');
    ylabel('PSF Attenuation Scale - LogRMS','FontSize',fontSize,'FontWeight','bold');
    legend(legendCellA,'Location','southeast'); legend boxoff; % columnlegend(3,legendCellA,'FontSize',24);
    set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on; axis square;
%     saveFigure(gcf, gca, [0,0,20,20], sprintf('Logarithmic_RMS_Gaussian_Effect_GL%d',glvect(indexGL)), 'svg');
end

%%  Arranged LogRMS per Gaussian witdh per TC
for indexGL = 1:length(glvect)
    perGLvariation = [AA(indexGL,:,1); AA(indexGL,:,2); AA(indexGL,:,3); AA(indexGL,:,4); AA(indexGL,:,5); AA(indexGL,:,6)];
    figure('color', 'white');
    hold on; arrayfun(@(index) plot(tcvect,perGLvariation(index,:),plotSpec{index},'color',colorSet(index,:),'LineWidth',lineWidth), plotRangeB); hold off;
    title(sprintf('Effect of Gaussian illumination-GL:%d',glvect(indexGL)),'FontSize',fontSize,'FontWeight','bold');
    xlabel('Topological Charge','FontSize',fontSize,'FontWeight','bold');
    ylabel('PSF Attenuation Scale - LogRMS','FontSize',fontSize,'FontWeight','bold');
    legend(legendCellB,'Location','southeast'); % columnlegend(3,legendCellA,'FontSize',24);
    set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on; axis square; xlim([min(tcvect), max(tcvect)]);
    %         saveFigure(gcf, gca, [0,0,20,20], sprintf('Logarithmic_RMS_Gaussian_N_4%d',1), 'svg');
end


%% Plots -- Arranged PowerSupression per Gaussian witdh per Gray Level
for indexGL = 1:length(glvect)
    PSperPSGLvariation = [BB(indexGL,:,1); BB(indexGL,:,2); BB(indexGL,:,3); BB(indexGL,:,4); BB(indexGL,:,5); BB(indexGL,:,6)]';
    figure('color', 'white');
    hold on; arrayfun(@(index) plot(gaussianWidth,PSperPSGLvariation(index,:),plotSpec{index},'color',colorSet(index,:),'LineWidth',lineWidth), plotRangeA); hold off;
    title(sprintf('Effect of Gaussian illumination-GL:%d',glvect(indexGL)),'FontSize',fontSize,'FontWeight','bold');
    xlabel('Gaussian Width [mm]','FontSize',fontSize,'FontWeight','bold');
    ylabel('EEF at Airy radius','FontSize',fontSize,'FontWeight','bold');
    legend(legendCellA,'Location','southeast'); % columnlegend(3,legendCellA,'FontSize',24);
    set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on; axis square;
    %         saveFigure(gcf, gca, [0,0,20,20], sprintf('Logarithmic_RMS_Gaussian_N_4%d',1), 'svg');
end

%%  Arranged PowerSupression per Gaussian witdh per TC
for indexGL = 1:length(glvect)
    PSperPSGLvariation = [BB(indexGL,:,1); BB(indexGL,:,2); BB(indexGL,:,3); BB(indexGL,:,4); BB(indexGL,:,5); BB(indexGL,:,6)];
    figure('color', 'white');
    hold on; arrayfun(@(index) plot(tcvect,PSperPSGLvariation(index,:),plotSpec{index},'color',colorSet(index,:),'LineWidth',lineWidth), plotRangeB); hold off;
    title(sprintf('Effect of Gaussian illumination-GL:%d',glvect(indexGL)),'FontSize',fontSize,'FontWeight','bold');
    xlabel('Topological Charge','FontSize',fontSize,'FontWeight','bold');
    ylabel('EEF at Airy radius','FontSize',fontSize,'FontWeight','bold');
    legend(legendCellB,'Location','northwest'); legend boxoff; % columnlegend(3,legendCellA,'FontSize',24);
    set(gca,'FontSize',fontSize,'FontWeight','normal'); grid on; axis square; xlim([min(tcvect), max(tcvect)]);
    saveFigure(gcf, gca, [0,0,20,20], sprintf('PowerSupression_Gaussian%d',glvect(indexGL)), 'svg');
end