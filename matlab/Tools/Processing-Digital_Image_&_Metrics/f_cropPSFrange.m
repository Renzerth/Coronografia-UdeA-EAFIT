function [croppedData] = f_cropPSFrange(PSFdata,cropRange)
%% Crop an image from a reference centered range
croppedData = PSFdata(cropRange(2,:),cropRange(1,:));
end