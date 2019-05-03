function f_plotGradient(cartcoord,GradEnergy,normIntensity,fontSize,lineWidth)
%% Plot gradient of intensity and intensity
figure('color','white');
[hAxes,hLine1,hLine2] = plotyy(cartcoord,GradEnergy,cartcoord,normIntensity);
set(hLine1, 'Color','b','LineStyle','-','LineWidth',lineWidth);
set(hLine2, 'Color','r','LineStyle','-','LineWidth',lineWidth);
set(hAxes(1),'YColor','b'); set(hAxes(2),'YColor','r');
set(hAxes,'FontSize',fontSize,'FontWeight','bold');
axis(hAxes(1), 'square'); axis(hAxes(2), 'square'); box(hAxes(1), 'on');

title('Analysis of Gradient"s intensity'); grid on;
xlabel('Angular position [\lambda/D]');
ylabel(hAxes(1),'Gradient of the intensity','FontSize',fontSize,'FontWeight','bold');
ylabel(hAxes(2),'Intensity pattern','FontSize',fontSize,'FontWeight','bold'); grid on;
lgdHandler = legend({'Gradient of the Intensity', 'Intensity'}); lgdHandler.FontSize = fontSize;
end