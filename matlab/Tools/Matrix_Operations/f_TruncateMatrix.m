function [Atrunc] = f_TruncateMatrix(B,BmidX,BmidY,A,AmidX,AmidY)
% Truncates matrix A so that it has the same size as B. The truncation
% occurrs in the center of A, whose midpoints can be found with
% f_ComputeMatrixMidPoints.m
% Inputs: 2D complex or real matrices A and B and midpoints of A and B
% Output: truncated matrix A
% Note: A must be greater than B
if (size(A,1) > size(B,1)) || (size(A,2) > size(B,2))
    newY = AmidY-BmidY+1 : AmidY+BmidY; % +1 since indices start with 1
    newX = AmidX-BmidX+1 : AmidX+BmidX; % +1 since indices start with 1
    Atrunc = A(newY,newX);
else
    disp(['Warning: A must be strictly bigger than B. A will not be' ...
         ' truncated']);
    Atrunc = A;
end

