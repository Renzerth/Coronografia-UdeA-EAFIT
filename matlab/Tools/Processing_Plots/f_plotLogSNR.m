function f_plotLogSNR(cartcoord,SNR,titprof,xlab)
%% Program Settings
tit = 'Log SNR';
tit1 = strcat(tit,{' '},titprof);
%% Plot of the Logarithmic SNR
figure('color','white');
plot(cartcoord,SNR);
xlabel(xlab); ylabel('SNR')
title(tit1); grid on;
end