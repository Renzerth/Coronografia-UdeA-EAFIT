function [exitPupilLens, cutOffFreq] = computeLens(Fx, Fy, diameter, focalLength, lambda)
% Fx, Fy are Spectral meshgrid coordinates
%
% In a thin lens system all Principal Planes and Pupils lies on the same plane
% exitPupilDiameter = lensDiameter; % The exit pupil diameter depends on the stop diameter
% exitPupilDistance = focalLengthA; %exit pupil distance equals to image distance in a thin lens system will be the EFL

relativeApertureFN = focalLength/diameter;
cutOffFreq =1/(lambda*relativeApertureFN); % [Cycles/m]
coordScale = 2/cutOffFreq; % Exit Pupil frequential scaling 

freqRadii = (Fx*coordScale).^2 + (Fy*coordScale).^2;
aperture = abs(sqrt(Fx.^2 + Fy.^2)*coordScale) <=1;
lensFase = exp(-1i*freqRadii.^2);

exitPupilLens = aperture.*lensFase;
end