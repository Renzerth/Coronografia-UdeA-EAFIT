function [midX,midY] = f_ComputeMatrixMidPoints(A)
% Computes the midpoints of a 2D matrix
% Inputs: 2D matrix (real or complex)
% Outputs: mid points in X and Y
MidF = @(S,idx) floor((size(S,idx)+1)/2); % Takes into account even and odd
                                          % size
midX = MidF(A,2); % X midpoint
midY = MidF(A,1); % Y midpoint

end

