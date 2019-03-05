%% Test if the center of the mask is inside the truncation
%   if plotMask == 2
%     tol = 0.003; % Sometimes it it not fully 0
%     [centY, centX] = find(rSLM>=0 & rSLM<=tol); % [row,col]. Finds the center of 
%                    % the polar radius so that the truncation is made on it
%                    % &: bitwise operator
%     if isempty(centY) || isempty(centX) % Error
%       error(['The center of the mask is out of the boundaries of the ' ...
%              'image, please select a smaller value for "shiftCart" ' ...
%              '(values in [0,100])']);
%     end
%   end