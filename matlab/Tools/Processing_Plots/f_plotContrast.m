function f_plotContrast(cartcoord,RefMeasurement,Measurement,titprof)
%% Program Settings
tit = strcat('Raw contrast',{' '},titprof);
%% Plot of each intensity in log scale
figure('color','white');
semilogy(cartcoord,Measurement); hold on
semilogy(cartcoord,RefMeasurement); hold off;
title(tit); grid on;
xlabel('Angular position [\lambda/D]'); ylabel('Relative contrast of the radial intensities [logscale]')
legend({'Coronographic', 'Non-Coronographic'});
end