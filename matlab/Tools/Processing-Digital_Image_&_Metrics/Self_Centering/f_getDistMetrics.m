function [criteriaValue, relativeChange, currentProfile] = getDistMetrics(currentFrame,dataSize, mainDataCenter, referenceRadialProfile, dataRange)
[currentProfile] = f_getAverageRadialProfile(currentFrame, dataSize, mainDataCenter);
profileLocalChanges = abs(referenceRadialProfile - currentProfile);
relativeChange = profileLocalChanges./referenceRadialProfile;
criteriaValue = mean(relativeChange(1:dataRange));
end