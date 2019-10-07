function [latexCharCell] = fetchDataToTable(varargin)
% Example modded to pass data and default parameters to the latexTable
%% Example N: using a splitted numerical array in different variables as data input
% clear input;
fprintf('\n\nExample N: using an array as data input\n\n');

% numeric values you want to tabulate:
% this field has to be an array or a MATLAB table
% in this example we use an array

if nargin == 0
    error('No variables were given.');
else
    rowNumber = length(varargin);
    colNumber = length(varargin{1});
    localData = zeros(rowNumber,colNumber);
    for varIndex = 1:rowNumber
        localData(varIndex,:) = varargin{varIndex};
    end
end
input.data = localData;

% Optional fields:

% Set column labels (use empty string for no label):
input.tableColLabels = transpose(strtrim(cellstr(num2str((1:colNumber)', 'col%d'))));

% Set row labels (use empty string for no label):
input.tableRowLabels = transpose(strtrim(cellstr(num2str((1:rowNumber)', 'row%d'))));

% Switch transposing/pivoting your table:
input.transposeTable = 0;

% Determine whether input.dataFormat is applied column or row based:
input.dataFormatMode = 'column'; % use 'column' or 'row'. if not set 'colum' is used

% Formatting-string to set the precision of the table values:
% For using different formats in different rows use a cell array like
% {myFormatString1,numberOfValues1,myFormatString2,numberOfValues2, ... }
% where myFormatString_ are formatting-strings and numberOfValues_ are the
% number of table columns or rows that the preceding formatting-string applies.
% Please make sure the sum of numberOfValues_ matches the number of columns or
% rows in input.tableData!
%
input.dataFormat = {'%.3f',colNumber}; % three digits precision for first two columns, one digit for the last

% Define how NaN values in input.tableData should be printed in the LaTex table:
input.dataNanString = '-';

% Column alignment in Latex table ('l'=left-justified, 'c'=centered,'r'=right-justified):
input.tableColumnAlignment = 'c';

% Switch table borders on/off (borders are enabled by default):
input.tableBorders = 1;

% Uses booktabs basic formating rules ('1' = using booktabs, '0' = not using booktabs).
% Note that this option requires the booktabs package being available in your LaTex.
% Also, setting the booktabs option to '1' overwrites input.tableBorders if it exists.
input.booktabs = 1;

% LaTex table caption:
input.tableCaption = 'MyTableCaption';

% LaTex table label:
input.tableLabel = 'MyTableLabel';

% Switch to generate a complete LaTex document or just a table:
input.makeCompleteLatexDocument = 1;

% call latexTable:
latexCharCell = latexTable(input);
end