%% Ranges
oddRange = 1:2:10;
evenRange = 2:2:10;
RMSFormula = @(dataSet, dataRange) sqrt(mean((mean(dataSet(:,dataRange)) - mean(mean(dataSet(:,dataRange)))).^2));

%% ODD Charges - GL Average LRMS per TC
fprintf('ODD Range - LRMS mean:\n');
LRMSmeanOdd = mean(arrangedLogRMS(:,oddRange)); disp(LRMSmeanOdd);
fprintf('ODD Range - LRMS mean STD:\n');
LRMSmeanSTDodd = std(arrangedLogRMS(:,oddRange)); disp(LRMSmeanSTDodd);

fprintf('ODD Range -  LRMS mean mean:\n');
LRMSglobalMeanOdd = mean(mean(arrangedLogRMS(:,oddRange))); disp(LRMSglobalMeanOdd);
fprintf('ODD Range - RMS error of LRMS mean:\n');
LRMSglobalMeanRMSodd = RMSFormula(arrangedLogRMS,oddRange); disp (LRMSglobalMeanRMSodd);

fprintf('ODD Range - Porcentual Improvement mean:\n');
GLImprovementMeanOdd = mean(GLlmprovement(:,oddRange)); disp(GLImprovementMeanOdd);
fprintf('ODD Range - Porcentual Improvement mean STD:\n');
GLImprovementMeanSTDodd = std(GLlmprovement(:,oddRange)); disp(GLImprovementMeanSTDodd);

fprintf('ODD Range - Power Supression mean:\n');
powerSupressOdd = transpose(mean(powerSupr(oddRange,:),2)); disp(powerSupressOdd);
fprintf('ODD Range - Power Supression mean STD:\n');
powerSupressSTDodd = std(transpose(powerSupr(oddRange,:))); disp(powerSupressSTDodd);

fprintf('\n\r');
%% EVEN Charges - GL Average LRMS per TC
fprintf('EVEN Range - LRMS mean:\n');
LRMSmeanEven = mean(arrangedLogRMS(:,evenRange)); disp(LRMSmeanEven);
fprintf('EVEN Range - LRMS mean STD:\n');
LRMSmeanSTDeven = std(arrangedLogRMS(:,evenRange)); disp(LRMSmeanSTDeven);

fprintf('EVEN Range -  LRMS mean mean:\n');
LRMSglobalMeanEven = mean(mean(arrangedLogRMS(:,evenRange))); disp(LRMSglobalMeanEven);
fprintf('EVEN Range - RMS error of LRMS mean:\n');
LRMSglobalMeanRMSeven = RMSFormula(arrangedLogRMS,evenRange); disp (LRMSglobalMeanRMSeven);

fprintf('EVEN Range - Porcentual Improvement mean:\n');
GLImprovementMeanEven= mean(GLlmprovement(:,evenRange)); disp(GLImprovementMeanEven);
fprintf('EVEN Range - Porcentual Improvement mean STD:\n');
GLImprovementMeanSTDeven = std(GLlmprovement(:,evenRange)); disp(GLImprovementMeanSTDeven);

fprintf('EVEN Range - Power Supression mean:\n');
powerSupressEven = transpose(mean(powerSupr(evenRange,:),2)); disp(powerSupressEven);
fprintf('EVEN Range - Power Supression mean STD:\n');
powerSupressSTDeven = std(transpose(powerSupr(evenRange,:))); disp(powerSupressSTDeven);

%% Make Latex Table 
[latexCharCellOdd] = fetchDataToTable(LRMSmeanOdd,LRMSmeanSTDodd,GLImprovementMeanOdd,GLImprovementMeanSTDodd,powerSupressOdd,powerSupressSTDodd);
[latexCharCellEven] = fetchDataToTable(LRMSmeanEven,LRMSmeanSTDeven,GLImprovementMeanEven,GLImprovementMeanSTDeven,powerSupressEven,powerSupressSTDeven);