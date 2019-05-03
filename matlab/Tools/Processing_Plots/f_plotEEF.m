function f_plotEEF(cartcoord,energy,normIntensity,titprof,xlab,fontSize,lineWidth,colorSet,lineStyle,markerSet)
%% Program Settings
tit = 'Encircled Energy Distribution of Intensity';
%% Plot comparison Energy and EEF distributions
figure('color','white');
[hAxes,hLine1,hLine2] = plotyy(cartcoord,energy,cartcoord,normIntensity);
set(hLine1, 'Color',colorSet(1,:),'LineStyle',lineStyle,'LineWidth', ...
            lineWidth,'Marker',markerSet);
set(hLine2, 'Color',colorSet(2,:),'LineStyle',lineStyle,'LineWidth', ...
            lineWidth,'Marker',markerSet);
set(hAxes(1),'YColor',colorSet(1,:)); set(hAxes(2),'YColor',colorSet(2,:));
set(hAxes,'FontSize',fontSize,'FontWeight','bold');

axis(hAxes(1), 'square'); axis(hAxes(2), 'square'); box(hAxes(1), 'on');
title(strcat(tit,{' '},titprof),'FontSize',fontSize,'FontWeight','bold'); grid on;
xlabel(xlab,'FontSize',fontSize,'FontWeight','bold');
ylabel(hAxes(1),'Relative throughput','FontSize',fontSize,'FontWeight','bold');
ylabel(hAxes(2),'Intensity pattern','FontSize',fontSize,'FontWeight','bold'); grid on;
lgdHandler = legend({'Encircled Energy Factor', 'Intensity'}); lgdHandler.FontSize = fontSize;
end