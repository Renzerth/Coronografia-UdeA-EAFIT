function [x,y,Hprof,Vprof,maxX,maxY,midX,midY] = f_makeImageProfile(x,y, ...
           dataArray,tol,shiftCart,tit,xlab,ylab,plotData,plotH,plotV,oneSideProfile,dcShift)
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
%  xlab,ylab: x and y labels (strings)
%  plotData: plots dataArray and lines of the profiles
%  plotH,plotV: booleans to plot the horizontal and the vertical profiles
%  oneSideProfile: 0 (no); 1(1st and 4th quadrants); 2(2nd and 3rd quadrants)
%  dcShift: Modulus applied for Fourier mask since the fftshift moves one
%  pixel the dc component. Only used for spectra (Fourier analysis)
%
% Outputs:
%  [Hprof,Vprof]: array of the profiles

%% Max points of the array (with the shift application)
% Takes into account an even/odd size of the size of dataArray
[maxY, maxX] = size(dataArray); % Size of the array

%% Mid points of the array
midX = round((maxX+1)/2) + dcShift*mod(maxX,2); % x mid point
midY = round((maxY+1)/2) + dcShift*mod(maxY,2); % y mid point

%% Shift of the mask scaling
shiftX = shiftCart(2); % Cartesian shift in x
shiftY = shiftCart(1); % Cartesian shift in y
shiftX = round(midX*shiftX); % Percentage of the half support of the mask
shiftY = round(midY*shiftY); % Percentage of the half support of the mask

%% Mid points of the array (with the shift application)
% The extreme parts of the profile are conserved and the midpoints shift
% The signs here are the opposite of the Shift application section of the
% function "f_DefineSpace" since they compensate this shift
midX = midX + shiftX; % Shifted X for the SLM 
midY = midY - shiftY; % Shifted Y for the SLM

%% Profile coordinates
Hx = [1,maxX]; % Horizontal x components: (xi,xf)
Vy = [1,maxY]; % Vertical y components: (yi,yf)
Hy = [midY,midY]; % Horizontal y components: (yi,yf)
Vx = [midX,midX]; % Vertical x components: (xi,xf)

%% Profiles of the 2D array drawn on the dataArray 
if plotData
    figure; imagesc(x,y,dataArray); title(tit); % colormap(hot)
    xlabel(xlab);
    ylabel(ylab);
    hold on
    % Horizontal:
    line(x(Hx),y(Hy),'LineWidth',3,'Color','blue','LineStyle','--'); 
    % Vertical:
    line(x(Vx),y(Vy),'LineWidth',3,'Color','red','LineStyle','--'); 
    hold off
end

%% Calculation of the profiles
Hprof = improfile(dataArray,Hx,Hy); % Horizontal profile
Vprof = improfile(dataArray,Vx,Vy); % Vertical profile

%% One side profile extraction
switch oneSideProfile
    case 0
        % No change
    case 1 %  Profiles in the 1st and 4th quadrants
        Hprof = Hprof(midX:end); % From the middle to the end
        Vprof = Vprof(midY:end); % From the middle to the end
        x = x(midX:end); % From the middle to the end
        y = y(midY:end); % From the middle to the end
    case 2% Profiles in the 2nd and 3rd quadrants
        Hprof = Hprof(1:midX); % From the middle to the end
        Vprof = Vprof(1:midY); % From the middle to the end
        x = x(1:midX); % From the middle to the end
        y = y(1:midY); % From the middle to the end        
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