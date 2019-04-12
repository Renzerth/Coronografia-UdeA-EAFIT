%% Cartesian coordinates
close all
x= -1:2/255:1; y = -1:2/512:1;
halfx = length(x)/2;
halfy = length(y)/2;
[X,Y] = meshgrid(x,y);
f = @(xi,yi) exp(0.01*xi.^2+0.02*yi.^2);
%% Vortex Mask
[Xr,Yr,aspectRatio,monitorSize] = f_MakeScreenCoords(3,false);
pupilRadius = 0.25/aspectRatio;
scaledY = Yr/aspectRatio;
maskGen = @(X,Y,r0,TC) mat2gray((sqrt(X.^2 + Y.^2) <= r0).* ...
                 angle(exp(1i*TC*(atan2(Y,X) + pi/TC))));

%% Example figure
% F = f(X,Y);
F = maskGen(Xr,scaledY,pupilRadius,10);
a = figure; b = imagesc(F);

%% Manual slider adjustment
[shiftYcoord, shiftXcoord] = f_addSliderPositioning(a,b);

%% All quadrants
X = X-shiftXcoord/halfx;
Y = Y-shiftYcoord/halfy;

%% Test figure with the shifts
% G=f(X,Y);
G = maskGen(Xr-shiftXcoord/(monitorSize(1)/2),scaledY-shiftYcoord/(aspectRatio*(monitorSize(1)/2)),pupilRadius,10);
figure; imagesc(x,y,G);
% xlim([-1 1]); ylim([-1 1]);

%%
monitorSize = [1920, 1080];
halfMx = monitorSize(1)/2;
halfMy = monitorSize(2)/2;

manualShiftYcoord = 0;
manualShiftXcoord = 90;
localX = 0;
localY = 1620;

shiftXi = localX + manualShiftXcoord; % Referred Overall SHIFTS from screen origin
shiftYi = +(localY + manualShiftYcoord); % Circshift Y periodicity correction
if shiftYi <= halfMy % shiftYi is never negative
  coorY = shiftYi;
else
  coorY = sign(shiftYi-halfMy)*(shiftYi + monitorSize(2)*(shiftYi <= halfMy) - monitorSize(2)) - halfMy*(halfMy == shiftYi); % Circshift Y periodicity correction 
%   shiftY = shiftYi + monitorSize(2)*(abs(shiftYi) > halfMy);
%   coorY = shiftY - halfMy;
  
end

if shiftXi <= halfMx % shiftXi is never negative
  coorX = shiftXi;
else
  coorX = sign(shiftXi-halfMx)*(shiftXi + monitorSize(1)*(shiftXi <= halfMx) - monitorSize(1)) - halfMx*(halfMx == shiftXi); % Circshift X periodicity correction 
  % coorX = shiftX;
  % OLD coorX = shiftX - halfMx;
end

