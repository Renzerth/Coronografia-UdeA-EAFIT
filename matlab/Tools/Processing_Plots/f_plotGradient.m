function f_plotGradient(cartcoord,GradEnergy,normIntensity)
%% Plot gradient of intensity and intensity
figure('color','white');
[hAxes,hLine1,hLine2] = plotyy(cartcoord,GradEnergy,cartcoord,normIntensity);
set(hLine1, 'Color','b','LineStyle','-'); set(hLine2, 'Color','r','LineStyle','-');
set(hAxes(1),'YColor','b'); set(hAxes(2),'YColor','r'); set(hAxes,'FontSize',12,'FontWeight','bold');
axis(hAxes(1), 'square'); axis(hAxes(2), 'square'); box(hAxes(1), 'on');
title(strcat(tit,{' '},titprof)); grid on;
xlabel('Angular position [\lambda/D]');
ylabel(hAxes(1),'Gradient of the relative throughput','FontSize',12,'FontWeight','bold');
ylabel(hAxes(2),'Intensity pattern','FontSize',14,'FontWeight','bold'); grid on;
lgdHandler = legend({'Gradient of the Encircled Energy Factor', 'Intensity'}); lgdHandler.FontSize = 10;
end