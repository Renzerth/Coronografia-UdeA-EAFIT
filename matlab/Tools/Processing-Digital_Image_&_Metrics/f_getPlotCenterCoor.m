function [refPoints] = f_getPlotCenterCoor(dataSize,midX,midY,shiftCart)
%% Shift of the mask scaling
shiftX = shiftCart(2); % Cartesian shift in x
shiftY = shiftCart(1); % Cartesian shift in y

%% Mid points of the array (with the shift application)
midX = midX + shiftX; % Shifted X for the SLM
midY = midY - shiftY; % Shifted Y for the SLM

%% Profile coordinates
Hx = [1,dataSize(1)]; % Horizontal x components: (xi,xf)
Vy = [1,dataSize(2)]; % Vertical y components: (yi,yf)
Hy = [midY,midY]; % Horizontal y components: (yi,yf)
Vx = [midX,midX]; % Vertical x components: (xi,xf)

%% Reference Data
refPoints = [Hx;Vy;Hy;Vx];
end