%% path dependencies
addpath(genpath(fileparts(pwd)));

%% Filter properties
maxTC = 10;
% GLRanges = 2:10;
GLRanges = [12 16 24 32 64 128 256];

%% Propagate PSF & Lyot Distributions
[LyotPlaneIntensities, PSFplaneIntensities, LyotReference, PSFreference] = computePSFVariationsAS_mod(maxTC, GLRanges);

%% Process PSF Intensities
f_ProcessSimulatedPSF();

%% Process Lyot Intensities
f_ProcessSimulatedLyot();
