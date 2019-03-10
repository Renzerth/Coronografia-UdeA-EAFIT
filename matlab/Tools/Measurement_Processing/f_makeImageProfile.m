function [Hprof,Vprof] = f_makeImageProfile(x,y,dataArray,tol,tit, ...
                         plotData,plotH,plotV,oneSideProfile)
% returns the profiles of a 2D data array (or matrix), i.e., an image or a
% 2D map. These profiles are horizontal and vertical
% Inputs:
%  x,y: cartesian coordinates for dataArray
%  dataArray: 2D numerical array (or matrix): image/mask for example
%  tol: discrete index tolerance in order to symmetrically dicrease the
%  size of the profile
%  tit: title for dataArray
%  plotData: plots dataArray and lines of the profiles
%  plotH,plotV: booleans to plot the horizontal and the vertical profiles
%  oneSideProfile:
% Outputs:
% [Hprof,Vprof]: array of the profiles

%% Mid and max points of the array
% Takes into account an even/odd size of the size of dataArray
[maxX, maxY] = size(dataArray); % Size of the array
midX = round((maxX+1)/2) + mod(maxX,2); % x mid point
midY = round((maxY+1)/2) + mod(maxY,2); % y mid point
 
%% Profiles of the 2D array drawn on the dataArray 
if plotData
    figure; imagesc(dataArray); title(tit); % colormap(hot)
    hold on
    % Horizontal:
    line([1,maxY],[midX,midX],'LineWidth',3,'Color','blue','LineStyle','--'); 
    % Vertical:
    line([midY,midY],[1,maxX],'LineWidth',3,'Color','red','LineStyle','--'); 
    hold off
end

%% Calculation of the profiles
Hprof = improfile(dataArray,[1,maxY],[midX,midX]); % Horizontal profile
Vprof = improfile(dataArray,[midY,midY],[1,maxX]); % Vertical profile

%% Plot of each profiles
if plotH
    %% Horizontal profile
    figure;
    plot(x(1+tol:end-tol),Hprof(1+tol:end-tol)); 
    title(['Horizontal profile of the ' tit]); 
end
if plotV
    %% Vertical profile
    figure;
    plot(y(1+tol:end-tol),Vprof(1+tol:end-tol)); 
    title(['Vertical profile of the ' tit]);
end 

end