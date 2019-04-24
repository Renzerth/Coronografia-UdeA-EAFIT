%% Profile of the measurements
for idxgral = 1:totalImgs
  [measEEFcurves{idxgral}] = f_calculateEEF(radialIntensityMeas{idxgral},cartcoord,titprof,tit,xlab);
end
%%
totalGL = length(glvect);
totalTC = length(tcvect);
powerSupr = zeros(totalTC,totalGL);
arrangedData = cell(totalTC,totalGL);
for tcIndx = 1:totalTC
  for glIndx = 1:totalGL
    idxgral = glIndx + (tcIndx - 1)*totalGL; % Reversed width/index
    powerSupr(tcIndx,glIndx) = measEEFcurves{idxgral}(aproxRadius);
    arrangedData{tcIndx,glIndx} = measEEFcurves{idxgral};
  end
end
%%
figure('color', 'white');
hold on; arrayfun(@(index) plot(glvect,powerSupr(index,:)), 1:totalTC); hold off;
xlabel('Gray Level'); ylabel('EEF at Airy Range');
legendCell = cellstr(num2str((1:10)', 'TC=%d')); legend(legendCell); grid on; axis square;
%%
plotAlotFunc = @(reference, Data) plot(reference,Data);
for indexGL = 1:totalGL
  figure();
  hold on; arrayfun(@(indexTC)plotAlotFunc(cartcoord,arrangedData{indexTC,indexGL}),1:totalTC); hold off;
  grid on; axis square; xlabel('Radial Distance (\lambda/D)'); ylabel('Throughput (EEF)'); legend(legendCell);
  title(sprintf('Throughput of Topological Charges at GL = %d',glvect(indexGL)));
end
%%
% cellfun
%
% [SNR] = f_calculateSNR(radialIntensityMeas{idxgral},radialIntensityRef,cartcoord,tit1,tit2,xlab);
%     case 3 % MSE
%       tit = 'Mean Squared Error';
%
%     otherwise
%       error('Select a valid metric');
%   end