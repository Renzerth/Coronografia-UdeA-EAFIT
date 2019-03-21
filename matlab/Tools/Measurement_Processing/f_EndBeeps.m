function f_EndBeeps(N,beepSound)
% N: number of beeps 
% beepSound: boolean for the sound
 if beepSound == 1
     for beepTimes = 1:N % Numbe of beeps
         beep();
         pause(0.2); % Time between the beeps
     end
 end
end

