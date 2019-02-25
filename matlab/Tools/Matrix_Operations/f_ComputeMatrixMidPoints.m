function [midX,midY] = f_ComputeMatrixMidPoints(A)
% Computes the midpoints of a 2D matrix
% Inputs: 2D matrix (real or complex)
% Outputs: mid points in X and Y
MidP = @(S,idx) floor((size(S,idx)+1)/2); % Takes into account even and odd
                                          % size
                                          % Should be ceiling
midX = MidP(A,2); % X midpoint
midY = MidP(A,1); % Y midpoint

end