function [shiftYcoord, shiftXcoord] = f_addSliderPositioning(parentHandler, figureHandler, varargin)

%% Input checking
if ~ishandle(parentHandler)
    error('Invalid or deleted figure handler.');
end

%% Program Settings
if nargin == 4 && ishandle(varargin{2})
  outsiderHandle = varargin{2};
elseif nargin == 3 && ishandle(varargin{1})
  outsiderHandle = varargin{1};
  sliderSize = 0.025; 
elseif nargin == 3 && isnumeric(varargin{1})
    sliderSize = varargin{1};
    outsiderHandle = [];
else
    outsiderHandle = [];
    sliderSize = 0.025; 
end
%% Slider Properties
dataXsize = figureHandler.XData(2);
dataYsize = figureHandler.YData(2);
stepRangeX = 1./[dataXsize,dataXsize];
stepRangeY = 1./[dataXsize,dataXsize];

%% Reference Points
halfPointX = (dataXsize+1)/2;
halfPointY = (dataYsize+1)/2;
originalData = figureHandler.CData;

%% Uicontrol properties
halfSlidSize = sliderSize/2;
set(parentHandler,'doublebuffer','on');

%% Uicontrol creation
xSliderShift = uicontrol('parent',parentHandler,'style','slider','units','normalized','position',[0,0,1,sliderSize],'min',0,'max',dataXsize,'value',halfPointX,'SliderStep',stepRangeX);
ySliderShift = uicontrol('parent',parentHandler,'style','slider','units','normalized','position',[1-halfSlidSize,sliderSize,halfSlidSize,1-sliderSize],'min',0,'max',dataYsize,'value',halfPointY,'SliderStep',stepRangeY);

%% User Slider operation
if ishandle(outsiderHandle)
  conditionHandle = outsiderHandle;
else
  conditionHandle = parentHandler;
end
while ishandle(conditionHandle) && ishandle(parentHandler) % Avoid usage of
                                        % addlistener to read slider values
    shiftYcoord = floor(ySliderShift.Value - halfPointY);
    shiftXcoord = floor(xSliderShift.Value - halfPointX);
    set(figureHandler,'CData',circshift(originalData,floor([shiftYcoord, shiftXcoord])));
    pause(0.1); % Allow idle instead of processing
    % Delay for user response: dead time of the processor
end


%% Former Ideas
% MyCallBack = @(a,b) disp(b.AffectedObject.Value);
% MyCallBack2 = @(src,event) disp(event.Key);
end