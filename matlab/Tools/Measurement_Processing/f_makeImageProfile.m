function [Hprof,Vprof,maxX,maxY,midX,midY] = f_makeImageProfile(x,y, ...
           dataArray,tol,shiftCart,tit,plotData,plotH,plotV,oneSideProfile)
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
%
% Outputs:
%  [Hprof,Vprof]: array of the profiles

%% Max points of the array (with the shift application)
% Takes into account an even/odd size of the size of dataArray
[maxX, maxY] = size(dataArray); % Size of the array


%% Mid points of the array
midX = round((maxX+1)/2); %+ mod(maxX,2); % x mid point
midY = round((maxY+1)/2); %+ mod(maxY,2); % y mid point


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

if shiftY >= 0 % To avoid it to be out of bounds
    maxY = maxY - shiftY; % Shifted Y for the SLM  
end

%% Profile coordinates
% The plus one takes into account that we deal with indices

if shiftX >= 0
    Hi = [1,maxX]; % Initial horizontal (x,y) pair
else
    Hi = [shiftX+1,maxX]; % Initial horizontal (x,y) pair
end

if shiftY >= 0
    Vf = [1,maxY]; % Final vertical (x,y) pair
else
    Vf = [shiftY+1,maxY]; % Final vertical (x,y) pair
end

Hf = [midY,midY]; % Final horizontal (x,y) pair
Vi = [midX,midX]; % Initial vertical (x,y) pair

%% Profiles of the 2D array drawn on the dataArray 
if plotData
    figure; imagesc(x,y,dataArray); title(tit); % colormap(hot)
    hold on
    % Horizontal:
    line(x(Hi),x(Hf),'LineWidth',3,'Color','blue','LineStyle','--'); 
    % Vertical:
    line(y(Vi),y(Vf),'LineWidth',3,'Color','red','LineStyle','--'); 
    hold off
end

%% Calculation of the profiles
Hprof = improfile(dataArray,Hi,Hf); % Horizontal profile
Vprof = improfile(dataArray,Vi,Vf); % Vertical profile

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