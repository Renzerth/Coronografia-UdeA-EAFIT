%% Ranges
oddRange = 1:2:10;
evenRange = 2:2:10;

%% ODD Charges - GL Average LRMS per TC 
mean(arrangedLogRMS(:,oddRange))
std(arrangedLogRMS(:,oddRange))

mean(mean(arrangedLogRMS(:,oddRange)))
sqrt(mean((mean(arrangedLogRMS(:,oddRange)) - mean(mean(arrangedLogRMS(:,oddRange)))).^2))

%% EVEN Charges - GL Average LRMS per TC 
mean(arrangedLogRMS(:,evenRange))
std(arrangedLogRMS(:,evenRange))

mean(mean(arrangedLogRMS(:,evenRange)))
sqrt(mean((mean(arrangedLogRMS(:,evenRange)) - mean(mean(arrangedLogRMS(:,evenRange)))).^2))