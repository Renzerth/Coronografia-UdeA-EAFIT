function [coorX, coorY] = f_circShiftToCart(circShiftX,circShiftY,shiftXPeriod,shiftYPeriod)
%% Coordinate Reference - Center
halfMx = shiftXPeriod/2;
halfMy = shiftYPeriod/2;
%% Correct Periodicity and refer to central cartesian coordinates
if circShiftY <= halfMy % shiftYi is never negative
  coorY = circShiftY;
else
  coorY = sign(circShiftY-halfMy)*(circShiftY + shiftYPeriod*(circShiftY <= halfMy) - shiftYPeriod) - halfMy*(halfMy == circShiftY); % Circshift Y periodicity correction 
end

if circShiftX <= halfMx % shiftXi is never negative
  coorX = circShiftX;
else
  coorX = sign(circShiftX-halfMx)*(circShiftX + shiftXPeriod*(circShiftX <= halfMx) - shiftXPeriod) - halfMx*(halfMx == circShiftX); % Circshift X periodicity correction 
end
end
