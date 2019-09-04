%% addpath dependencies
addpath(genpath(fileparts(pwd)));
%% Filter properties
maxTC = 10;
GLRanges = 2:10;
% GLRanges = [12 16 24 32 64 128 256];
%% Propagate PSF & Lyot Distributions
[LyotPlaneIntensities, PSFplaneIntensities, LyotReference, PSFreference] = computePSFVariationsAS(maxTC, GLRanges);
save(strcat(fileparts(pwd),'/Data/AngularSpectrumSimulations/angular_spectrum_pipeline.mat'))
%% Process PSF Intensities
% f_ProcessSimulatedPSF();
%% Process Lyot Intensities
f_ProcessSimulatedLyot();