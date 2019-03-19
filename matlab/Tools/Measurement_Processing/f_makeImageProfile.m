function [Hprof,Vprof,maxX,maxY,midX,midY] = f_makeImageProfile(x,y, ...
           dataArray,tol,shiftCart,tit,plotData,plotH,plotV,oneSideProfile,dcShift)
% returns the profiles of a 2D data array (or matrix), i.e., an image or a
% 2D map. These profiles are horizontal and vertical
%
% Inputs:
%  x,y: cartesian coordinates for dataArray
%  dataArray: 2D numerical array (or matrix): image/mask for example
%  tol: discrete index tolerance in order to symmetrically dicrease the
%  size of the profile
%  shiftCart:
%  tit: title for dataArray
%  plotData: plots dataArray and lines of the profiles
%  plotH,plotV: booleans to plot the horizontal and the vertical profiles
%  oneSideProfile:
%  dcShift: Modulus applied for Fourier mask since the fftshift moves one
% pixel the dc component. Only used for spectra (Fourier analysis)
%
% Outputs:
%  [Hprof,Vprof]: array of the profiles

%% Max points of the array (with the shift application)
% Takes into account an even/odd size of the size of dataArray
[maxX, maxY] = size(dataArray); % Size of the array


%% Mid points of the array
midX = round((maxX+1)/2) + dcShift*mod(maxX,2); % x mid point
midY = round((maxY+1)/2) + dcShift*mod(maxY,2); % y mid point


%% Shift of the mask scaling
shiftX = shiftCart(2); % Cartesian shift in x
shiftY = shiftCart(1); % Cartesian shift in y
shiftX = round(midX*shiftX); % Percentage of the half support of the mask
shiftY = round(midY*shiftY); % Percentage of the half support of the mask

%% Mid points of the array (with the shift application)
midX = midX + shiftX; % Shifted X for the SLM 
midY = midY - shiftY; % Shifted Y for the SLM

%% Max points of the array (with the shift application)
if shiftX <= 0 % To avoid it to be out of bounds
    maxX = maxX + shiftX; % Shifted X for the SLM 
end

if shiftY <= 0 % To avoid it to be out of bounds
    maxY = maxY + shiftY; % Shifted Y for the SLM  
end

%% Profile coordinates
% The plus one takes into account that we deal with indices
% if shiftX <= 0
    Hx = [1,maxX]; % Horizontal x components: (xi,xf)
% else
%     Hx = [shiftX+1,maxX]; % Horizontal x components: (xi,xf)
% end

% if shiftY <= 0
    Vy = [1,maxY]; % Vertical y components: (yi,yf)
% else
%     Vy = [shiftY+1,maxY]; % Vertical y components: (yi,yf)
% end

Hy = [midY,midY]; % Horizontal y components: (yi,yf)
Vx = [midX,midX]; % Vertical x components: (xi,xf)

%% Profiles of the 2D array drawn on the dataArray 
if plotData
    figure; imagesc(x,y,dataArray); title(tit); % colormap(hot)
    hold on
    % Horizontal:
    line(x(Hx),x(Hy),'LineWidth',3,'Color','blue','LineStyle','--'); 
    % Vertical:
    line(y(Vx),y(Vy),'LineWidth',3,'Color','red','LineStyle','--'); 
    hold off
end

%% Calculation of the profiles
Hprof = improfile(dataArray,Hx,Hy); % Horizontal profile
Vprof = improfile(dataArray,Vx,Vy); % Vertical profile

%% One side profile extraction
if oneSideProfile
    Hprof = Hprof(fix(end/2):end); % From the middle to the end
    Vprof = Vprof(fix(end/2):end); % From the middle to the end
    x = x(fix(end/2):end); % From the middle to the end
    y = y(fix(end/2):end); % From the middle to the end
end 

%% Plot of each profiles
if plotH
    %% Horizontal profile
    figure;
    plot(x(1+tol:end-tol),Hprof(1+tol:end-tol)); 
    title(strcat('Horizontal profile of the',{' '},tit)); 
end
if plotV
    %% Vertical profile
    figure;
    plot(y(1+tol:end-tol),Vprof(1+tol:end-tol)); 
    title(strcat('Vertical profile of the',{' '},tit));
end 

end