%% Ranges
oddRange = 1:2:10;
evenRange = 2:2:10;

%% ODD Charges - GL Average LRMS per TC 
mean(arrangedLogRMS(:,oddRange))
std(arrangedLogRMS(:,oddRange))

mean(mean(arrangedLogRMS(:,oddRange)))
sqrt(mean((mean(arrangedLogRMS(:,oddRange)) - mean(mean(arrangedLogRMS(:,oddRange)))).^2))

mean(GLlmprovement(:,oddRange))
std(GLlmprovement(:,oddRange))

transpose(mean(powerSupr(plotRange,:),2))
transpose(std(powerSupr(plotRange,:),2))

%% EVEN Charges - GL Average LRMS per TC 
mean(arrangedLogRMS(:,evenRange))
std(arrangedLogRMS(:,evenRange))

mean(mean(arrangedLogRMS(:,evenRange)))
sqrt(mean((mean(arrangedLogRMS(:,evenRange)) - mean(mean(arrangedLogRMS(:,evenRange)))).^2))

mean(GLlmprovement(:,evenRange))
std(GLlmprovement(:,evenRange))

transpose(mean(powerSupr(evenRange,:),2))
transpose(std(powerSupr(evenRange,:),2))

%%