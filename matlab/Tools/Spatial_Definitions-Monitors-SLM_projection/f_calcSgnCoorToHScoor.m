function [SgnWidth, SgnHeight] = f_calcSgnCoorToHScoor(width, height)
% This function does the opposite of the following:
% (oposite operation of "f_calcHScoorToSgnCoor")
% changes origin refereded coordinates width = [0,W] and height = [0,H]
% that are normalized in ranges of [0,1] to a proportional range [-1,1].
% Input coordinates can be treated as a Heavyside function with a domain
% of [0,1]. To change input range between [-width/2,width/2] and 
% [-height,height] in a proportional range [-1,1] the relationship between
% the Sign function and the Heavyside function is used.
%
% H = 0.5(1+Sgn(x)) :: Domain [0,1] <- Sgn(x) = 2H - 1 :: Domain [-1,1]

%% Check input
if abs(width) > 1 || abs(height) > 1
  error('Height and Width Values must be defined between [-1,1]');
end

%% Transform Coordinates
SgnWidth = 0.5*(1 + width);
SgnHeight = 0.5*(1 + height);

end