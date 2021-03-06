function Apadded = f_PadMatrix(A,B,method)
% Pads the matrix A in order to fit the matrix B by adding values to A
% A must be smaller than B
  if (size(B,1) > size(A,1)) || (size(B,2) > size(A,2))  
    c = size(B) - size(A); 
    idx = find(c ~= 0);
    pad = c(idx);
    hor = 1; vert = 2; % Padding Direction: horizontal (H) or vertical (V)
    if isscalar(pad)                  
     % Replicated since that the outer value of the Zernike are constant
     if idx == 2 % Horizontal Padding 
     Apadded = f_MatrixSymmetricPadding(A,pad,method,hor); % Horizontal
     else % idx = 1, Vertical Padding
      Apadded = f_MatrixSymmetricPadding(A,pad,method,vert); % Vertical
     end
    else % pad is a vector
      Apadded = f_MatrixSymmetricPadding(A,pad(1),method,vert); % Vertical
      Apadded = f_MatrixSymmetricPadding(Apadded,pad(2),method,hor); % Hor
    end
  else
    disp(['Warning: A must be strictly smaller than B. A will not' ... 
         ' be padded']);
    Apadded = A;
  end
end

