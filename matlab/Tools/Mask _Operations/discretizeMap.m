function [phi] = discretizeMap(maskMap,NG)
phaseVor = mod(maskMap, 2*pi);  
phi = 2*pi*floor(phaseVor/(2*pi/NG))/NG;
end