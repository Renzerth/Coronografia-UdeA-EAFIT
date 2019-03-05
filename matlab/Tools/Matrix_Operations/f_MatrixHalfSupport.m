function HalfSupportA = f_MatrixHalfSupport( A )
% Computes the half support of a matrix A that is supposed to be a
% symmetric meshgrid

absminA = abs(min( A(:) ));
absmaxA = abs(max( A(:) ));
HalfSupportA = (absminA + absmaxA)/2;
end

