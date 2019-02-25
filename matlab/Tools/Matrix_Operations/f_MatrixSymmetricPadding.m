function Apadded = f_MatrixSymmetricPadding(A,pad,method,direction)
% Calculates a symmetric padarray and takes into account wether pad is even
% or odd
% 
% Inputs:
%  A: matrix to be padded
%  pad: size of the padding (scalar). It is divided between both sides of
%       the matrix. Must be bigger than 0
%  method: 'replicate', 'symmetric', 'circular' or a scalar
%  direction: (1): horizontal; (2): vertical
%
% Output:
%  Apadded: padded matrix

if pad >= 0
  midP = @(n) ceil(n/2); % mid point function

  % 'Horizontal' case:
  a1 = [0 midP(pad)]; % (both, for even pad) OR (pre, for odd pad)
  a2 = [0 midP(pad) - 1]; % pos, for odd pad only
  if direction == 2 % 'Vertical' case:
      a1 = circshift(a1,[0,1]); % Circular shift
      a2 = circshift(a2,[0,1]); % Circular shift
  end

  %% Even or odd padding
  if mod(pad,2) == 0 % Even pad (symmetric)
      Apadded = padarray(A,a1,method,'both');
  else % Odd pad (asymmetric)
      Apadded = padarray(A,a1,method,'pre');
      Apadded = padarray(Apadded,a2,method,'pos');
  end  
else
  Apadded = A;
  warning('pad must be bigger than zero');
end
end