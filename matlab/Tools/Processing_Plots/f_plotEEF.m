function f_plotEEF(cartcoord, energy, normIntensity, titprof, xlab)
%% Program Settings
tit = 'Encircled Energy Distribution of Intensity';
%% Plot comparison Energy and EEF distributions
figure('color','white');
[hAxes,hLine1,hLine2] = plotyy(cartcoord,energy,cartcoord,normIntensity);
set(hLine1, 'Color','b','LineStyle','-'); set(hLine2, 'Color','r','LineStyle','-');
set(hAxes(1),'YColor','b'); set(hAxes(2),'YColor','r'); set(hAxes,'FontSize',10,'FontWeight','bold');
axis(hAxes(1), 'square'); axis(hAxes(2), 'square'); box(hAxes(1), 'on');
title(strcat(tit,{' '},titprof)); grid on;
xlabel(xlab);
ylabel(hAxes(1),'Relative throughput','FontSize',14,'FontWeight','bold');
ylabel(hAxes(2),'Intensity pattern','FontSize',14,'FontWeight','bold'); grid on;
lgdHandler = legend({'Encircled Energy Factor', 'Intensity'}); lgdHandler.FontSize = 10;
end