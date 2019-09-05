%% path dependencies
addpath(genpath(fileparts(pwd)));

%% Filter properties
maxTC = 10;
% GLRanges = 2:10;
GLRanges = [12 16 24 32 64 128 256];

%% Propagate PSF & Lyot Distributions
[LyotPlaneIntensities, PSFplaneIntensities, LyotReference, PSFreference] = computePSFVariationsAS(maxTC, GLRanges);
try
    save(strcat(fileparts(pwd),'/Data/AngularSpectrumSimulations/angular_spectrum_pipeline.mat'))
catch
    TGTdir = strcat(fileparts(pwd),'/Data/AngularSpectrumSimulations');
    mkdir(TGTdir)
    save(strcat(TGTdir,'/angular_spectrum_pipeline.mat'))
end

%% Process PSF Intensities
f_ProcessSimulatedPSF();

%% Process Lyot Intensities
f_ProcessSimulatedLyot();
