function [criteriaValue,relativeChange,currentProfile] = ...
f_getDistMetrics(currentFrame,dataSize, mainDataCenter, ...
referenceRadialProfile,dataRange)
[currentProfile] = f_getAverageRadialProfile(currentFrame,dataSize, ...
                                             mainDataCenter);
profileLocalChanges = abs(referenceRadialProfile - currentProfile);
relativeChange = profileLocalChanges./referenceRadialProfile;
criteriaValue = mean(relativeChange(dataRange));

end