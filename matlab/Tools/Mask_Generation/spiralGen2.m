function [vortexMask] = spiralGen2(fieldSize,TC)
% Centered vortex mask generation

m = fieldSize(1); % Row Size
n = fieldSize(2); % Column Size

Row = (1:m) -1-m/2; % All Rows
Column1 = (1:(n/2)+1) - 1 -n/2; % First left half
Column2 = (n/2+1:n) - 1 -n/2; % Last right half

[x1,y1] = meshgrid(Column1,Row);
th1 = atan(y1./x1) + pi/2; % First Half Angular transition

[x2,y2] = meshgrid(Column2,Row);
th2 = atan(y2./x2) + pi+pi/2; % Last Half Angular transition

angularTransition = [th1(:,1:(n/2)),th2(:,1:end)];

vortexMask = exp(1i*TC*angularTransition);
vortexMask(m/2+1,n/2+1) = 0; % Remove atan(0/0) nan value
end