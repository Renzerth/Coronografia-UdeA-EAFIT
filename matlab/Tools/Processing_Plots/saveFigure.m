function saveFigure(figureHandler, figureAxes, figurePosition, fileName, fileFormat)
%% Define Figure Size
set(figureHandler, 'Units', 'centimeters','Position',figurePosition)
%% Tight Axes
tightWhiteSpace(figureAxes)
%% Set Paper Properties
figureHandler.PaperPositionMode = 'auto';
figPosition = figureHandler.PaperPosition;
figureHandler.PaperSize = [figPosition(3) figPosition(4)];
%% Save figure
print(figureHandler,fileName,sprintf('-d%s',fileFormat));
end