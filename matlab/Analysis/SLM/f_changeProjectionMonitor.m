function [targetResolution,monitorsInformation] = f_changeProjectionMonitor(varargin)
%% Input checking
if nargin == 1 && ~isnumeric(varargin{1})
    if ismember(varargin{1},{'Restore','restore','R','r'})
        set(0,'DefaultFigurePosition','default');
        targetResolution = [];
        monitorsInformation = [];
        return;
    else
        error('Monitor index must be an integer value.');
    end
elseif nargin == 1 && isnumeric(varargin{1})
    screenIndex = varargin{1};
    selectMaxRes = false;
    enableChange = true;
elseif nargin == 2 && ~isnumeric(varargin{1})
    error('Monitor index must be an integer value.');
elseif nargin == 2 && isnumeric(varargin{1})
    screenIndex = varargin{1};
    selectMaxRes = false;
    enableChange = varargin{2};
else
    selectMaxRes = true;
    enableChange = true;
end
%% Monitor properties
monitorsInformation = get(0,'MonitorPositions');
primaryMonitorInfo = get(0,'ScreenSize');
availableMonitors = monitorsInformation(~ismember(monitorsInformation,primaryMonitorInfo,'rows'),:);
availableResolutions = availableMonitors(:,3:end);
locationShifts = availableMonitors(:,1:end/2);
%% Verify available monitor
if size(monitorsInformation,1) < screenIndex
    warning('Selected monitor is not available. Maximum resolution will be selected.');
    selectMaxRes = true;
end
%% Define target monitor
if selectMaxRes
    [targetResolution,screenIndex] = max(availableResolutions,[],1);
    targetShifts = locationShifts(screenIndex(1),:);
    if enableChange
        set(0,'DefaultFigurePosition',[targetShifts targetResolution]);
    end
else
    targetResolution = monitorsInformation(screenIndex,3:end);
    if enableChange
        set(0,'DefaultFigurePosition',monitorsInformation(screenIndex,:));
    end
end

end