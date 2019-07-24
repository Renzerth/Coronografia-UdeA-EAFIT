function f_plotMosaic(dataArrange,xRefVector,yRefVector,titleSet,xLabelSet,yLabelSet,customMap,fontSize,saveEnabled,abs_ang,enableAxis)
% Patrick Martineau
% Perfect subplot in Matlab
% http://p-martineau.com/perfect-subplot-in-matlab/
%% Initialization
[subplotsx,subplotsy] = size(dataArrange);
%% Parameters for figure and panel size
plotheight = 20; % cm
plotwidth = 16;

leftedge = 1.2;
rightedge = 0.4;
topedge = 1;
bottomedge = 1.5;

spacex = 0.2;
spacey = 0.2;
fontsizeAx = fontSize - 5; % Axes font size
fontSizeTit = fontSize - 3; % Labels & Titles font size

%% Compute placement properties
sub_pos = subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

%% Setting the Matlab figure
f=figure('visible','on','color','white');
clf(f);
set(gcf, 'PaperUnits', 'centimeters');
set(gcf, 'PaperSize', [plotwidth plotheight]);
set(gcf, 'PaperPositionMode', 'manual');
set(gcf, 'PaperPosition', [0 0 plotwidth plotheight]);

%% Loop to create axes
for i=1:subplotsx
    for ii=1:subplotsy
        
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsizeAx,'Box','on','Layer','top','FontWeight','bold');
        imagesc(xRefVector,yRefVector,dataArrange{i,ii}); colormap(customMap);
        cbarHandler=colorbar; limVals=get(cbarHandler,'Limits');
        %         if abs_ang == 1
        tol2 = 0.2*limVals(2); % Colorbar custom tick adjustment
        set(cbarHandler,'Ticks',linspace(limVals(1)+tol2,limVals(2)-tol2,2));
        % NO: ticksLabels = cellstr(num2str(limVals', '%1.0e')); SCIENTIFIC
        ticksLabels = cellstr(num2str(limVals', '%1.3f'));
        set(cbarHandler,'XTickLabel',ticksLabels);
        %         else
        %             yticks([-0.8 0 0.8]); xticks([-0.8 0 0.8]);
        %             xticklabels({'-1','0','1'});yticklabels({'-1','0','1'});
        %             tol2 = 0.1*pi; % Colorbar custom tick adjustment
        %             caxis([-pi,pi]);
        %             set(cbarHandler,'Ticks',linspace(-pi+tol2,pi-tol2,2));
        %             set(cbarHandler,'TickLabels',{'-\pi' '\pi'}); % For the masks
        %         end
        
        if ii==subplotsy
            title(titleSet{i},'FontSize',fontSizeTit,'FontWeight','bold');
        end
        
        if ii>1
            set(ax,'xticklabel',[]);
        end
        
        if i>1
            set(ax,'yticklabel',[]);
        end
        
        if i==1
            ylabel(yLabelSet{ii},'FontSize',fontSizeTit,'FontWeight','bold');
        end
        
        if ii==1
            xlabel(xLabelSet{i},'FontSize',fontSizeTit,'FontWeight','bold');
        end
        set(gca,'FontSize',fontSize-1);
    end
end

%% Square axis
axesHandles = findobj(get(gcf,'Children'), 'flat','Type','axes');
axis(axesHandles,'square');

if isa(enableAxis,'logical') && ~enableAxis
    axis(axesHandles,'off')
end

%% Saving eps with matlab and then producing pdf and png with system commands
if saveEnabled
    filename='test';
    print(gcf, '-depsc2','-loose',[filename,'.eps']);
    system(['epstopdf ',filename,'.eps'])
    system(['convert -density 300 ',filename,'.eps ',filename,'.png'])
end
end