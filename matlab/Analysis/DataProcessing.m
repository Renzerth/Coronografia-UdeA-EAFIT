%% Post-processing of the data (applycation of the metric of the degree of extintion)
for idxgral = 1:totalImgs
    cd(DatalogDir); % Goes to the data log directory
    A = load(MeasInfo{idxgral},'A');
    cd(ProcessedDir); % Goes to the output directory
    A = A'; % Processing of the image
    save(['processed_' MeasInfo{idxgral}],'A');
end