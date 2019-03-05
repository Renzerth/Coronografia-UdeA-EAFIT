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

  %% Radius for the circular (or elliptical) pupil truncation
  % In both cases of the next if-else, rMax is found as the maximum radius
  % that allows to circumscribe a circle inside a square or inside a 
  % rectangle
  
  % if coordType == 1 % User-defined
  %     rMax = max(r(:));  % the maximum value of r (diagonal of the square)
  %     rSize = rMax/sqrt(2); % Equals this since twice rSize^2 equals
  %                           % rmax^2 (Pythagorean theorem)
  % else % coordType == 2 % Screen-resolution defined
  % end

%   rSizeVect = size(r); % 2D vector with the size of r
%   rMidVect = floor((rSizeVect+1)/2); % mid points of the size of r
%   rMove = max(rMidVect); % Maximum midpoint in order to move here
%   idx = find(rMidVect == rMove); % finds the index or rMove to determine
%                                  % if one should move in the x or the y 
%                                  % direction
%   if isscalar(idx)
%     if idx == 1 % (for landscape monitors)
%           rMax = r(rMove,1); % The rMax is in the y direction 
%     else % idx == 2 % (for portrait monitors)
%           rMax = r(1,rMove); % The rMax is in the x direction
%     end
%   else % idx is a two-row vector, meanning that one has square-sized figures
%     rMax = r(1,rMove)/sqrt(2);
%   end
%   rSize = rMax; % Both rmax and rsize are equal
                  % if one has a unitary space, rMax = 1 always
                  
%   rSizeVect = size(r);
%   rminSize = min(rSizeVect);
%   rmaxSize = max(rSizeVect);
%   rSizeIdx =  floor((rminSize+1)/2);
%   xtrunc = linspace(-1,1,rminSize);
%   ytrunc = linspace(-1,1,rmaxSize);
%   [Xtrunc,Ytrunc] = meshgrid(xtrunc,ytrunc);
%   [~,rtrunc] = cart2pol(Xtrunc,Ytrunc);
%   rSize = min(rtrunc(rSizeIdx,1),rtrunc(1,rSizeIdx)); 
  % Meaning: min(landscapeMonitor,portraitMonitor): the selected minimum
  % takes into account both possible screen configurations
%   figure; imagesc(r);
  
 % rSize = 1; % Works for CoordType = 2