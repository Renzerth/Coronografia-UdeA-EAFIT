function f_plotAttenRatios(cartcoord,attenuationData,leakageData,fontSize,lineWidth)
%% Plot Attenuation Ratios of Intensity
figure('color','white');
[hAxes,hLine1,hLine2] = plotyy(cartcoord,attenuationData,cartcoord,leakageData);
set(hLine1, 'Color','b','LineStyle','-','LineWidth',lineWidth); 
set(hLine2, 'Color','r','LineStyle','-','LineWidth',lineWidth);
set(hAxes(1),'YColor','b'); set(hAxes(2),'YColor','r');
set(hAxes,'FontSize',fontSize,'FontWeight','bold');
axis(hAxes(1), 'square'); axis(hAxes(2), 'square'); box(hAxes(1), 'on');

title('Photonic Attenuation Ratio of Intensity'); grid on;
xlabel('Angular position [\lambda/D]');
ylabel(hAxes(1),'Attenuation','FontSize',fontSize,'FontWeight','bold');
ylabel(hAxes(2),'Leakage','FontSize',fontSize,'FontWeight','bold'); grid on;
lgdHandler = legend({'Attenuation', 'Leakage'}); lgdHandler.FontSize = fontSize;
end