function [VZk,Pupil]=zernikes(p,apperture,m_size)
% zernikes is a function that computes the Zernike
% expansion base polynomials. This program features
% p-th polynomial calculation using the Noll index 
% notation. Also matrix size can be configured. 
% Normalization factor is Pi. But it is forcibly
% normalized to one to obtain an orthonormal
% expansion base
%
% Inputs:  p - Positive integer - Number of polynomials.
%             apperture - float - Defines pupil size
%             m_size - Positive integer - Matrix size.
%
% Outputs: VZk -  Vector of Zernike's Polynomials
%                       Used for Coefficient decomposition.
%              Pupil - NaN/Integer matrix - Unitary circle.
%
% Version: 1.2 for Matlab
%
% Author: Juan Jos? Cadavid Mu?oz. EAFIT University
%
% Date:14/03/2015.
%% Intializing
Pupil = NaN(m_size);
Nnm = NaN(1,p);
%% Index convertion
j = (0:p)'; % Array with polynomials to be used in the representation
n = ceil((-3+sqrt(9+8*j))/2); % Polynomial Order array
m = 2*j-n.*(n+2); % Azimuthal Frequency array
%% Integration element
Delta = (2/m_size/apperture)^2;
%% Space construction
[X,Y] = meshgrid(-1:2/(m_size-1):1);
[Phi, Rho] = cart2pol(X,Y);
%% Pupil creation
Pupil(Rho<apperture) = 1;
%% Unitary circle mask 
Pupil(Pupil==0) = NaN;
Pupil(Pupil==1) = 0;
%% Circle coordinates scaling
Rho=Rho/apperture;
%% Normalized Coordinates within unit circle
Rho=Rho(~isnan(Pupil));
Phi=Phi(~isnan(Pupil));
%% Store vector initializing
VZk = NaN(size(Pupil(~isnan(Pupil)),1),p);
nanPupilSize = size(Pupil(~isnan(Pupil)));
%% Zernike Polynomials calculation
for i=1:p
    Nnm(i)=sqrt(2*(n(i)+1)/(1+(m(i)==0))); % Normalization constant
    Radial = zeros(nanPupilSize);
    for s=0:0.5*(n(i)-abs(m(i)))    % Radial part calculation
        Radial = Radial + ((-1)^s.*factorial(n(i)-s).*(Rho).^(n(i)-2.*s)/...
        (factorial(s).*factorial(0.5*(n(i)+abs(m(i)))-s).*factorial(0.5.*(n(i)-abs(m(i)))-s)));
    end
    Angular= ((m(i)>=0).*cos(m(i).*Phi)-(m(i)<0).*sin(m(i).*Phi)); %Angular part calculation
    Radial=(Nnm(i).*Radial.*Angular);
    VZk(:,i) = Radial(:);
end
%% Unit normalization
Iprod = VZk.'* VZk*Delta;
VZk = bsxfun(@rdivide, VZk, sqrt(diag(Iprod).'));
%% Kronecker Delta Plotting
% zProdsNorm = VZk.' * VZk * Delta;
% figure, imagesc(zProdsNorm), colorbar;
end