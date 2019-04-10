manualShiftYcoord = 0;
manualShiftXcoord = 0;
localX = shiftXcoord;
localY = shiftYcoord;

monitorSize = [1920, 1080];

shiftX = localX + manualShiftXcoord; % Referred Overall SHIFTS from screen origin
shiftX = sign(shiftX-monitorSize(1)/2)*(shiftX + monitorSize(1)*(shiftX <= monitorSize(1)/2) - monitorSize(1)) - monitorSize(1)/2*(monitorSize(1)/2==shiftX); % Circshift X periodicity correction
shiftY = -(localY + manualShiftYcoord); % Circshift Y periodicity correction
shiftY = shiftY + monitorSize(2)*(abs(shiftY) > monitorSize(2)/2);

coorX = monitorSize(1)/2 + shiftX;
coorY = shiftY - monitorSize(2)/2;
