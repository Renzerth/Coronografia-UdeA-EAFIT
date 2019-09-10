function [angularTransition] = compAngTransition(fieldSize,varargin)
% Centered vortex mask generation

if nargin == 2 && isnumeric(varargin{1}) && length(varargin{1}) == 2
    shiftY = varargin{1}(1); % % Apply center shift when available
    shiftX = varargin{1}(2);
else
    shiftY = 0; % By default no shift is made and center is half size +1
    shiftX = 0;
end

m = fieldSize(1); % Row Size
n = fieldSize(2); % Column Size

halfRow = floor((m+1)/2);
halfCol = floor((n+1)/2);

centerPointRow = (halfRow + shiftY); % displace row center
centerPointCol = (halfCol + shiftX); % displace col center

Row = (1:m) - (centerPointRow+1); % All Rows
Column1 = (1:(centerPointCol)+1) - (centerPointCol+1); % First left half
Column2 = (centerPointCol+1:n) - (centerPointCol+1); % Last right half

[x1,y1] = meshgrid(Column1,Row);
th1 = atan(y1./x1) + pi/2; % First Half Angular transition

[x2,y2] = meshgrid(Column2,Row);
th2 = atan(y2./x2) + pi+pi/2; % Last Half Angular transition

angularTransition = [th1(:,1:(centerPointCol)),th2(:,1:end)];
angularTransition(centerPointRow+1,centerPointCol+1) = 2*pi;

end