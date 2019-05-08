function f_plotLogSNR(cartcoord,SNR,titprof,xlab,fontSize,lineWidth)
%% Program Settings
tit = 'Log SNR';
tit1 = strcat(tit,{' '},titprof);
%% Plot of the Logarithmic SNR
figure('color','white');
plot(cartcoord,SNR,'LineWidth',lineWidth);
xlabel(xlab,'FontSize',fontSize,'FontWeight','bold');
ylabel('SNR','FontSize',fontSize,'FontWeight','bold');
title(tit1,'FontSize',fontSize,'FontWeight','bold'); grid on;
set(gca,'FontSize',fontSize,'FontWeight','normal');
end