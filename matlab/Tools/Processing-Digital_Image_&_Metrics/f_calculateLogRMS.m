function [logRMS] = f_calculateLogRMS(logSNR)
%% Logarithmic Coronographic - Non-coronographic deviation
logSNR(isinf(logSNR)) = nan; % Drop inf values of zero-valued intensities
logRMS = sqrt(nanmean(logSNR.^2));
end