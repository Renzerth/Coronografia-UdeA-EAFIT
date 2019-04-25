function f_plotSNR(cartcoord,RefMeasurement,Measurement,SNR,titprof,xlab)
%% Program Settings
tit = 'Signal-to-Noise Ratio';
tit1 = strcat(tit,{' '},titprof);
tit2 = strcat('Raw contrast',{' '},titprof);
%% Plot of each intensity in log scale
figure('color','white');
semilogy(cartcoord,Measurement); hold on
semilogy(cartcoord,RefMeasurement); hold off;
title(tit2); grid on;
xlabel('Angular position [\lambda/D]'); ylabel('Relative contrast of the radial intensities [logscale]')
legend({'Measurement', 'Reference'});

%% Plot of the Logarithmic SNR
figure('color','white');
plot(cartcoord,SNR);
xlabel(xlab); ylabel('SNR')
title(tit1); grid on;
end