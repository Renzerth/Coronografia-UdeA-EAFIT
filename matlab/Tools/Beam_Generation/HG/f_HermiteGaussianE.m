%----------------------------------------------------------------------------------------------
% PROGRAM: HermiteGaussianE
% AUTHOR:  Andri M. Gretarsson
% DATE:    6/26/03
%
% SYNTAX: z=HermiteGaussianE([l,m,q <,lambda,a>],x <,y>);
%           <...> indicates optional arguments
%
% Returns the complex field amplitude of a single polarization Hermite-Gaussian mode.
% The domain of the function is specified in cartesian coordinates x,y. 
% If the input argument a is not specified or is equal to 1, the area under the
% surface z(x,y) is equal to 1.  If z=z(x), i.e. the one dimensional form is being used,
% then the area under the surface formed by rotating z(x) around the z axis equals 1.
% In short the functions are normalized and form an orthonormal basis wrt. the inner
% product formed by multiplying two HermiteE's together and integrating over area.
% The formula is taken from "Principles of Lasers" by Orazio Svelto, 4th ed. Section 5.5.1.3.
%
% Note that it is possible to specify the parameters l,m,q,a,lambda as column vectors.  In that case,
% the matrix returned z, will be a stack of 2D matrixes (i.e. a 3D matrix), where successive 2D matrixes
% in the stack (i.e. successive layers of the 3D matrix) correpsond to successive rows of the 
% l,m,q,a,lambda columns
%
% INPUT ARGUMENTS:
% ----------------
% l,m    = Hermite Gaussian mode numbers (integers, positive or zero). Can be a scalar or a column vector.
% q      = complex radius of curvature of the beam (1/q = 1/R + i*lambda/pi/w^2). 
%          Can be a scalar or a column vector.
% a      = complex prefactor ( includes phase, e.g. for a beam that
%          has been propageted with an ABCD matrix a = 1/(A + [B/q1])^(1+l+m) ).
%          See for example, O. Svelto, Principles of Lasers 4th ed. eqn. 4.7.30.
%          Can be a scalar or a column vector.
% lambda = Wavelength of the light. Can be a scalar or a column vector.
% x      = x position vector or a 2D matrix generated by meshgrid. 
% y      = y position vector or a 2D matrix generated by meshgrid.
%          If y is not specified then the default y is used:  
%          if x is 1D, y=zeros(size(x)). If x is 2D, y=x.
%
% OUTPUT ARGUMENTS:
% -----------------
% z(i,j) = Complex field of the Hermite Gaussian mode with parameters l,m,... at 
%          ( x(i,j),y(i,j) ).  May be a vector or a matrix depending on whether 
%          x and y are vectors or matrixes.  If in addition, l,m,... are columns, 
%          then z(i,j,k) will correspond to the complex field amplitude of the 
%          Hermite Gaussian mode corresponding to the parameters l(k),m(k),...
%
% NOTES:
% ------
% If x and y are not vectors but matrixes generated by
% the matlab function meshgrid, then the output variable z
% is a matrix rather than a vector.  The matrix form allows
% the function HermiteGaussianE to have a plane as it's domain 
% rather than a curve.
%
% If the parameters l,m,q,p,lambda are equal length column vectors rather 
% than scalars, z is a size(x)*length(lambda) matrix.  E.g. if size(x) is n*n
% then each level z(:,:,k) is a 2D field of a Hermite Gaussian with
% the parameters given by [l(k),m(k),q(k),lambda(k),a(k)].
%
% EXAMPLE 1 (2D):  
%      w=[0.001; 0.001];           
%      xseed=[-3*max(w):max(w)/30:3*max(w)];  yseed=xseed';
%      [x,y]=meshgrid(xseed,yseed);
%      lambda = [1.064e-6 ; 1.064e-6];
%      R = [-30 ; -30];
%      q = (1./R - i* lambda./pi./w.^2).^(-1);  a=[1;1];
%      l=[1;2]; m=[0;2];
%      E=HermiteGaussianE([l,m,q,lambda,a],x,y);
%      subplot(2,1,1); h1=pcolor(x,y,abs(E(:,:,1)));  set(h1,'EdgeColor','none'); axis square;
%      subplot(2,1,2); h2=pcolor(x,y,abs(E(:,:,2)));  set(h2,'EdgeColor','none'); axis square; shg;
%
% EXAMPLE 2 (1D):
%       w=[1,2,3,4].'; x=[-10:0.001:10]; lambda=[1,1,1,1].'*656e-9; R=[1,1,1,1].'*1000; 
%       q = (1./R - i* lambda./pi./w.^2).^(-1);  a=[1,1,1,1].'; l = [0,1,2,3].'; m=[0,0,0,0].';
%       E=HermiteGaussianE([l,m,q,lambda,a],x); I=E.*conj(E); phi=angle(E);
%       plot(x,I(:,1),x,I(:,2),x,I(:,3),x,I(:,4));
%
% Last updated: July 18, 2004 by AMG
%----------------------------------------------------------------------------------------------
%% SYNTAX: z = f_HermiteGaussianE([l,m,q <,lambda,a>],x <,y>);
%----------------------------------------------------------------------------------------------

function z = f_HermiteGaussianE(params,x,varargin)

defaultcoord2=0;
if nargin>=3, y=varargin{1}; else defaultcoord2=1;  end

if min(size(x))==1
    if size(x,1)<size(x,2), x=x'; end  %make x and y columnar
    if defaultcoord2
        y=zeros(size(x));
    else 
        if size(y,1)<size(y,2), y=y'; end
    end
        
end

if min(size(x)) > 1 
    z=zeros(size(x,1),size(x,2),size(params,1));  % need this since zeros(size(y),10) gives a 2D matrix even if y is 2D!  (Matlab feature.)
    if defaultcoord2, y=transpose(x); end
else
    z=zeros(size(x),size(params,1)); 
end

l=params(:,1);
m=params(:,2);
q=params(:,3);
if size(params,2)>=4
    lambda=params(:,4);
else
    lambda=1064e-9;
end
if size(params,2)>=5
    a=params(:,5);
else
    a=ones(size(q));
end

w = f_axial_w_(q,lambda);

if min(size(x))>=2
    for u=1:size(params,1)
            z(:,:,u) = a(u)...
                .* sqrt((2/pi))/2^((l(u)+m(u))/2)/sqrt(factorial(l(u))*factorial(m(u)))/w(u)...
                .* f_HermitePoly(l(u),sqrt(2).*x/w(u)).*f_HermitePoly(m(u),sqrt(2).*y/w(u))...
                .* exp( -1i*2*pi/lambda(u)*(x.^2+y.^2)/2/q(u) ) ;
    end
else
    for u=1:size(params,1)
            z(:,u) = a(u)...
                .* sqrt((2/pi))/2^((l(u)+m(u))/2)/sqrt(factorial(l(u))*factorial(m(u)))/w(u)...
                .* transverse.f_HermitePoly(l(u),sqrt(2).*x/w(u)).*transverse.f_HermitePoly(m(u),sqrt(2).*y/w(u))...
                .* exp( -1i*2*pi/lambda(u)*(x.^2+y.^2)/2/q(u) );
    end
end