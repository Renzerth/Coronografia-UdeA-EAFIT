function [ hFig ] = f_CustomPreview(vid,tit)
%   Detailed explanation goes here
% from: https://nl.mathworks.com/help/imaq/previewing-data.html
% and: https://nl.mathworks.com/help/vision/examples/video-display-in-a-custom-user-interface.html
% and: https://nl.mathworks.com/help/imaq/examples/working-with-properties.html
%
% Taken from those links and organized by Samuel Plazas

% Create a figure window. This example turns off the default
% toolbar, menubar, and figure numbering:
figname = char(strcat('Preview of the ',{' '},'camera'));
hFig = figure('Toolbar','figure','Menubar', 'none','NumberTitle', ...
  'Off','Name',figname);

% Set up the push buttons:
uicontrol('String', 'Close','Callback', 'close(gcf)','Units',...
  'normalized','Position',[0 0 0.15 .07]);

% Create the text label for the timestamp:
hTextLabel = uicontrol('style','text','String','Timestamp', ...
  'Units','normalized','Position',[0.85 -.04 .15 .08]);

% Create the image object in which you want to display the video
% preview data.
vidRes = vid.VideoResolution;
imWidth = vidRes(1);
imHeight = vidRes(2);
nBands = vid.NumberOfBands;
hImage = image( zeros(imHeight, imWidth, nBands) );

% Set video title using uicontrol. uicontrol is used so that text
% can be positioned in the context of the figure, not the axis.
titlePos = [0.24 0.93 0.5 0.06];
axisTitle = tit;
uicontrol('style','text',...
  'String', axisTitle,...
  'Units','Normalized',...
  'Parent',hFig,'Position', titlePos,...
  'BackgroundColor',hFig.Color);

% Specify the size of the axes that contains the image object
% so that it displays the image at the right resolution and
% centers it in the figure window.
figSize = get(hFig,'Position');
figWidth = figSize(3);
figHeight = figSize(4);
gca.unit = 'pixels';
gca.position = [((figWidth - imWidth)/2) ((figHeight - imHeight)/2) ...
  imWidth imHeight];

% Set up the update preview window function.
setappdata(hImage,'UpdatePreviewWindowFcn',@f_mypreview_fcn);

% Make handle to text label available to update function.
setappdata(hImage,'HandleToTimestampLabel',hTextLabel);

% Display the video data in your GUI:
preview(vid, hImage);
end

