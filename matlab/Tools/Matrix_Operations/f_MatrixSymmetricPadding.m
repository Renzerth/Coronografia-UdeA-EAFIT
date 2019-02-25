function Apadded = f_MatrixSymmetricPadding(A,pad,method,direction)
% Calculates a symmetric padarray and takes into account wether pad is even
% or odd
% 
% Inputs:
%  A: matrix to be padded
%  pad: size of the padding (scalar)
%  method: 'replicate', 'symmetric', 'circular' or a scalar
%  direction: (1): horizontal; (2): vertical
%
% Output:
%  Apadded: padded matrix

% 'Vertical' case:
a1 = [0 pad/2];
a2 = [0 (pad+1)/2];
a3 = [0 (pad+1)/2 - 1];
if direction == 1 % 'Horizontal' case:
    a1 = circshift(a1,1); % Circular shift
    a2 = circshift(a2,1); % Circular shift
    a3 = circshift(a3,1); % Circular shift
end

%% Even or odd padding
if mod(pad,2) == 0 % Even pad (symmetric)
    Apadded = padarray(A,a1,method,'both');
else % Odd pad (asymmetric)
    Apadded = padarray(A,a2,method,'pre');
    Apadded = padarray(Apadded,a3,method,'pos');
end  
end