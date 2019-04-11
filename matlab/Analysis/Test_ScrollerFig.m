%% Cartesian coordinates
close all
x= -1:2/255:1; y = -1:2/512:1;
halfx = length(x)/2;
halfy = length(y)/2;
[X,Y] = meshgrid(x,y);
f = @(xi,yi) exp(0.01*xi.^2+0.02*yi.^2);

%% Example figure
F=f(X,Y);
a = figure; b = imagesc(F);

%% Manual slider adjustment
[shiftYcoord, shiftXcoord] = f_addSliderPositioning(a,b);

%% All quadrants
X = X-shiftXcoord/halfx;
Y = Y-shiftYcoord/halfy;

%% Test figure with the shifts
G=f(X,Y);
figure; imagesc(x,y,G);
% xlim([-1 1]); ylim([-1 1]);

%%
manualShiftYcoord = 0;
manualShiftXcoord =  0;
localX = shiftXcoord;
localY = shiftYcoord;

monitorSize = [1920, 1080];

shiftX = localX + manualShiftXcoord; % Referred Overall SHIFTS from screen origin
shiftX = sign(shiftX-monitorSize(1)/2)*(shiftX + monitorSize(1)*(shiftX <= monitorSize(1)/2) - monitorSize(1)) - monitorSize(1)/2*(monitorSize(1)/2==shiftX); % Circshift X periodicity correction
shiftY = -(localY + manualShiftYcoord); % Circshift Y periodicity correction
shiftY = shiftY + monitorSize(2)*(abs(shiftY) > monitorSize(2)/2);

coorX = monitorSize(1)/2 + shiftX;
coorY = shiftY - monitorSize(2)/2;
