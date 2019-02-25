function [VZk,Pupil] = f_ZernikePolynomials(p,apperture,m_size)
% zernikes is a function that computes the Zernike
% expansion base polynomials. This program features
% p-th polynomial calculation using the Noll index 
% notation. Also, the matrix size can be configured. 
% Normalization factor is Pi. But it is forcibly
% normalized to one to obtain an orthonormal
% expansion base (and almost orthogonal since
% it is a discrete base)
%
% Inputs:  p - Positive integer - Number of polynomials.
%          apperture - float - Defines pupil relative size (w.r.t. m_size)
%          m_size - Positive integer - Matrix size.
%
% Outputs: VZk -  Vector of Zernike's Polynomials
%                 used for Coefficient decomposition.
%          Pupil - NaN/Integer matrix - Unitary circle.
%
% Version: 1.2 for Matlab
%
% Author: Juan Jos� Cadavid Mu�oz. EAFIT University
% Date: 14/03/2015
% Commented by Samuel Plazas Escudero on 2018/04/03

%% Intializing
Pupil = NaN(m_size); % m_size x m_size matrix
Nnm = NaN(1,p); % NaN: Not a Number

%% Index conversion
j = (0:p)'; % Array with polynomials to be used in the representation
n = ceil((-3+sqrt(9+8*j))/2); % Polynomial Order array
m = 2*j-n.*(n+2); % Azimuthal Frequency array

%% Integration element
Delta = (2/m_size/apperture)^2;

%% Space construction
[X,Y] = meshgrid(-1:2/(m_size-1):1); % Unitary space
[Phi, Rho] = cart2pol(X,Y); % Polar

%% Pupil creation
Pupil(Rho<apperture) = 1; % One inside a circle or radius apperture

%% Unitary circle mask 
Pupil(Pupil==0) = NaN;
Pupil(Pupil==1) = 0;

%% Circle coordinates scaling
Rho = Rho/apperture;

%% Normalized Coordinates within unit circle
Rho = Rho(~isnan(Pupil)); % Assigns what isn't NaN on Pupil
Phi = Phi(~isnan(Pupil));

%% Store vector initializing
VZk = NaN(size(Pupil(~isnan(Pupil)),1),p);

%% Zernike Polynomials calculation
% with line 28 (Index conversion), n(i)-abs(m(i)) is automatically not zero
% for even results
for i = 1:p % p polynomials are generated
    Nnm(i)= sqrt(2*(n(i)+1)/(1+(m(i)==0))); % Normalization constant
                              % m(i)==0: one for i = 0 and zero for i =~ 0
    Radial = zeros(size(Pupil(~isnan(Pupil))));
    for s = 0:0.5*(n(i)-abs(m(i)))    % Radial part calculation
        Radial = Radial + ((-1)^s.*factorial(n(i)-s).*(Rho).^(n(i)-2.*s)/...
        (factorial(s).*factorial(0.5*(n(i)+abs(m(i)))-s).*...
        factorial(0.5.*(n(i)-abs(m(i)))-s)));
    end 
    % Angular part calculation
    Angular = ((m(i)>=0).*cos(m(i).*Phi)-(m(i)<0).*sin(m(i).*Phi)); 
    % cos for positives and zero, -sin for negatives
    Radial = (Nnm(i).*Radial.*Angular); % Normalized radial (with Nnm)
    VZk(:,i) = Radial(:); % Stores the i-th polynomial
end

%% Unit normalization
Iprod = VZk.'* VZk*Delta; % Product of the polynomials, this is a "discrete
                          % integral", or better said, a dot product of the
                          % polynomials
% A powerful function that avoids using a for:
VZk = bsxfun(@rdivide, VZk, sqrt(diag(Iprod).')); 
% Diagonal is divided by its element. Its elements initially approximate to
% pi (and not equal to pi), since one has a discrete base and the 
% polynomials aren't entirely orthogonal

%% Kronecker Delta Plotting 
% Ones on diagonal and zero outside: means they are correctly normalized 
% (although they aren't fully orthogonal)
% zProdsNorm = VZk.' * VZk * Delta;
% figure, imagesc(zProdsNorm), colorbar;
end