function [ alpha ] = f_LambdaDToarcsec(alpha)
alpha = rad2deg(alpha); % Angular distance in degrees
alpha = alpha*3600; % alpha[arcsecond] = alpha[rad]*3600[arcsecond]/1[rad]
end

