function [rejection, attenuation] = f_calculateAttenuat(nonCoronagraphic,coronagraphic)
%% Compute proportion of with-vortex and vortex intensities
rejection = nonCoronagraphic./coronagraphic;
attenuation = 1/rejection;
end