function [phi] = discretizeMap(maskMap,NG)
phaseVor = mod(maskMap + pi, 2*pi);  
phi = floor(phaseVor/(2*pi/NG))/NG;
phi = exp(1i*(2*pi*phi - pi));
end