function f_plotContrast(cartcoord,RefMeasurement,Measurement,titprof,fontSize,lineWidth)
%% Program Settings
tit = strcat('Raw contrast',{' '},titprof);
%% Plot of each intensity in log scale
figure('color','white');
semilogy(cartcoord,Measurement,'LineWidth',lineWidth); hold on
semilogy(cartcoord,RefMeasurement,'LineWidth',lineWidth); hold off;
title(tit,'FontSize',fontSize,'FontWeight','bold'); grid on;
xlabel('Angular separation [\lambda/D]','FontSize',fontSize,'FontWeight','bold');
ylabel('Relative contrast of the radial intensities [logscale]','FontSize',fontSize,'FontWeight','bold')
legend({'Coronagraphic', 'Non-Coronagraphic'},'FontSize',fontSize,'FontWeight','bold');
set(gca,'FontSize',fontSize,'FontWeight','normal');
end