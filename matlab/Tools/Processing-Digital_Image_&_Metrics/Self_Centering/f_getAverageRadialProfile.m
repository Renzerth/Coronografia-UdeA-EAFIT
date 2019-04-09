function [averageRadialProfile] = f_getAverageRadialProfile(grayImage, ...
                                       dataSize,coordinatesCenter,varargin)
% Thanks to Image Analyst for the original Demo
% www.mathworks.com/matlabcentral/answers/276298-how-to-plot-the-radial-profile-of-a-2d-image
%% Program Settings
if nargin == 4 && islogical(varargin{1})
    plotsEnabled = varargin{1};
else
    plotsEnabled = false;
end
xCenter = coordinatesCenter(1);
yCenter = coordinatesCenter(2);
rows = dataSize(1);
columns = dataSize(2); 
%% Find out what the max distance will be by computing the distance to each corner
distanceToUL = sqrt((1-yCenter)^2 + (1-xCenter)^2);
distanceToUR = sqrt((1-yCenter)^2 + (columns-xCenter)^2);
distanceToLL = sqrt((rows-yCenter)^2 + (1-xCenter)^2);
distanceToLR= sqrt((rows-yCenter)^2 + (columns-xCenter)^2);
maxDistance = ceil(max([distanceToUL, distanceToUR, distanceToLL, distanceToLR]));

%% Allocate an array for the profile
profileSums = zeros(1, maxDistance);
profileCounts = zeros(1, maxDistance);

%% Scan the original image getting gray level, and scan edtImage getting distance
% Then add those values to the profile.
for column = 1 : columns
    for row = 1 : rows
        thisDistance = round(sqrt((row - yCenter)^2 + (column - xCenter)^2));
        if thisDistance <= 0
            continue;
        end
        profileSums(thisDistance) = profileSums(thisDistance) + double(grayImage(row, column));
        profileCounts(thisDistance) = profileCounts(thisDistance) + 1;
    end
end
% Divide the sums by the counts at each distance to get the average profile
averageRadialProfile = profileSums ./ profileCounts;

%% Plots
if plotsEnabled
    fontSize = 20;
    
    subplot(2, 2, 1);
    imshow(grayImage, []);
    axis on;
    
    title('Original Image', 'FontSize', fontSize);
    % Set up figure properties:
    % Enlarge figure to full screen.
    set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
    % Get rid of tool bar and pulldown menus that are along top of figure.
    set(gcf, 'Toolbar', 'none', 'Menu', 'none');
    % Give a name to the title bar.
    set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off')
    drawnow;
    
    subplot(2, 2, 2);
    imshow(grayImage, []);
    axis on;
    title('With circles along radius gridlines from plot below', 'FontSize', fontSize);
    hold on;
    line([xCenter, xCenter], [1, rows], 'Color', 'r', 'LineWidth', 2);
    line([1, columns], [yCenter, yCenter], 'Color', 'r', 'LineWidth', 2);
    
    % Let's compute and display the histogram.  Just for fun - it's not necessary though.
    [pixelCount, grayLevels] = imhist(grayImage);
    subplot(2, 2, 3);
    bar(grayLevels, pixelCount);
    grid on;
    title('Histogram of Original Image', 'FontSize', fontSize);
    xlim([0 grayLevels(end)]); % Scale x axis manually.
    drawnow; % Force it to paint the screen.
    
    subplot(2, 2, 4);
    plot(1:length(averageRadialProfile), averageRadialProfile, 'b-', 'LineWidth', 3);
    grid on;
    title('Average Radial Profile', 'FontSize', fontSize);
    xlabel('Distance from center', 'FontSize', fontSize);
    ylabel('Average Gray Level', 'FontSize', fontSize);
    
    % We want to have the circles over the image be at the same distances as the grid lines along the x axis.
    % Get the tick marks along the x axis
    ax = gca;
    xTickNumbers = ax.XTick;
    xTickNumbers(xTickNumbers == 0) = [];  % Zero is probably in there, so get rid of it.
    radii = xTickNumbers;
    centers = repmat([xCenter, yCenter], length(radii), 1);
    subplot(2, 2, 2);  % Switch to upper right image.
    hold on;
    viscircles(centers,radii); % Plot the circles over the image.
end