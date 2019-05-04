function f_plotMosaic(dataArrange,xRefVector,yRefVector,titleSet,xLabelSet,yLabelSet,customMap,fontSize,saveEnabled)
% Patrick Martineau
% Perfect subplot in Matlab
% http://p-martineau.com/perfect-subplot-in-matlab/
%% Initialization
[subplotsx,subplotsy] = size(dataArrange);
%% Parameters for figure and panel size
plotheight=20;
plotwidth=16;

leftedge=1.2;
rightedge=0.4;
topedge=1;
bottomedge=1.5;

spacex=0.2;
spacey=0.2;
fontsize=fontSize-5; % Axes font size
fontSize = fontSize - 3; % Labels & Titles font size

%% Compute placement properties
sub_pos=subplot_pos(plotwidth,plotheight,leftedge,rightedge,bottomedge,topedge,subplotsx,subplotsy,spacex,spacey);

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
        
        ax=axes('position',sub_pos{i,ii},'XGrid','off','XMinorGrid','off','FontSize',fontsize,'Box','on','Layer','top','FontWeight','bold');
        
        imagesc(xRefVector,yRefVector,dataArrange{i,ii}); colormap(customMap);
        
        if ii==subplotsy
            title(titleSet{i},'FontSize',fontSize,'FontWeight','bold');
        end
        
        if ii>1
            set(ax,'xticklabel',[]);
        end
        
        if i>1
            set(ax,'yticklabel',[]);
        end
        
        if i==1
            ylabel(yLabelSet{ii});
        end
        
        if ii==1
            xlabel(xLabelSet{i});
        end
        set(gca,'FontSize',fontsize,'FontWeight','bold');
    end
end

%% Saving eps with matlab and then producing pdf and png with system commands
if saveEnabled
    filename='test';
    print(gcf, '-depsc2','-loose',[filename,'.eps']);
    system(['epstopdf ',filename,'.eps'])
    system(['convert -density 300 ',filename,'.eps ',filename,'.png'])
end
end