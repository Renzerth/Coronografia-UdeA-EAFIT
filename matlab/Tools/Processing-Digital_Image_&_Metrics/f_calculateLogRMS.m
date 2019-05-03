function [logRMS] = f_calculateLogRMS(logSNR)
%% Logarithmic Coronographic Non-coronographic deviation
logRMS = sqrt(mean(logSNR.^2));
end