function [PB,X,Y]=f_Parabolic_Nondiffracting_Beam(L,N,parity,kt,a)
%{
Calculate a parabolic nondiffracting beam 

INPUT:
    L = transverse physical size of the X-Y space [-L,~L,-L,~L]
    N = number of sampling points 
    parity = parity of the beam; 0 = EVEN  C,  1 = ODD 
    kt = transverse wavevector
    a  = "parabolicity" of the beam [-Inf,Inf]
          
OUTPUT:
    PB = Parabolic Nondiffracting Beam
    X and Y = space matrices

EXAMPLE:   
        [PB,X,Y]=Parabolic_Nondiffracting_Beam(6*10^-3,501,0,10^4,3); figure;
        surf(X,Y,abs(PB)); shading interp; lighting phong; view(2); axis equal tight; xlabel('x','FontName','Helvetica','Fontsize',15); ylabel('y','FontName','Helvetica','Fontsize',15); ht=title('Parabolic Nondiffracting Beam','Fontsize',12);
        set(ht,'FontName','Helvetica','FontSize',12,'FontWeight','bold'); set(gca,'FontName','Helvetica','FontSize',10,'FontWeight','bold'); colormap(hot);
    
For more information about parabolic nondiffracting beam :
	
    "Parabolic nondiffracting optical wavefields," Miguel A. Bandres, J. C. Gutiérrez-Vega, and S. Chávez-Cerda
    Optics Letters, 29(1), 44-46 (2004) ( http://goo.gl/KhmqQY )

    "Observation of Parabolic nondiffracting optical fields," Carlos López-Mariscal, Miguel A. Bandres, S. Chávez-Cerda, and J. C. Gutiérrez-Vega
    Optics Express, 13(7), 2364-2369 (2005) ( http://goo.gl/rFhTMD )

AUTHOR: Miguel A. Bandres
               www.mabandres.com
%}

%% Variables
g=sqrt(2*kt);
[xi,eta,X,Y]=cart2parabolic(L,N);   % Constructed parabolic and Cartesian coordinate meshes 
% G=(1+parity)/(pi*sqrt(2))*abs(mfun('gamma',parity*2/4+1/4+1/2*1i*a))^2; 
G=(1+parity)/(pi*sqrt(2))*abs(gamma(parity*2/4+1/4+1/2*1i*a))^2; 

%% Parabolic Nondiffracting Beam
if parity
        P_eta = yO_parabolic(-a,g*eta);
        P_xi   = yO_parabolic(a,g*xi);
    else
        P_eta = yE_parabolic(-a,g*eta);
        P_xi    = yE_parabolic(a,g*xi);    
end

PB=G*P_eta.*P_xi;   % Parabolic Nondiffractive Beam

return;

%% Subfunctions
function [v,u,X,Y]=cart2parabolic(d,N) % Cartesian Coordinates to Parabolic Coordinates

even=0; if mod(N,2)==0;  N=N+1;  even=1; end;

[X,Y]=meshgrid(linspace(-d,d,N),linspace(-d,d,N));
Y=flipud(Y);

u=zeros(size(X)); v=zeros(size(X)); 

% Calculate First quadrant
    u(1:(N+1)/2,(N+1)/2:N)=sqrt(sqrt(X(1:(N+1)/2,(N+1)/2:N).^2+Y(1:(N+1)/2,(N+1)/2:N).^2)+X(1:(N+1)/2,(N+1)/2:N));
    v(1:(N+1)/2,(N+1)/2:N)=sqrt(sqrt(X(1:(N+1)/2,(N+1)/2:N).^2+Y(1:(N+1)/2,(N+1)/2:N).^2)-X(1:(N+1)/2,(N+1)/2:N));

% Complete quadrants by symmetry 
    u(1:(N+1)/2,1:(N+1)/2)=fliplr(v(1:(N+1)/2,(N+1)/2:N));
    v(1:(N+1)/2,1:(N+1)/2)=fliplr(u(1:(N+1)/2,(N+1)/2:N));
    u((N+1)/2:N,1:N)=flipud(-u(1:(N+1)/2,1:N));
    v((N+1)/2:N,1:N)=flipud(v(1:(N+1)/2,1:N));

if even==1  
    u(N,:)=[]; u(:,N)=[];
    v(N,:)=[]; v(:,N)=[];
     X(N,:)=[]; X(:,N)=[];
     Y(N,:)=[]; Y(:,N)=[];        
end;

return;

function [y,yp]=yE_parabolic(a,z)       % EVEN Weber functions (parabolic cylinder functions)   
%{
This function calculates the EVEN Weber functions (parabolic cylinder functions) 
which are solutions of the Weber equation (parabolic cylinder equation):

d^2y/dz^2+(z^2/4-a)y=0

where
  -Infinity<z<Infinity 
   a=real parameter

INPUTS:
-Infinity<z<Infinity  REAL independent variable (Vector or Matrix). 
a = real parameter

OUTPUTS:
y =  EVEN Weber function
yo = First derivative of the EVEN Weber function

For more information about Weber functions:

    Miguel A. Bandres and B.M. Rodriguez-Lara, "Nondiffracting accelerating waves: Weber waves and parabolic momentum"
    New Journal of Physics, 15(013054) (2013) ( http://goo.gl/XfaVqq )

    Miguel A. Bandres, J. C. Gutiérrez-Vega, and S. Chávez-Cerda, "Parabolic nondiffracting optical wavefields"
    Optics Letters, 29(1), 44-46 (2004) ( http://goo.gl/KhmqQY )

AUTHOR: Miguel A. Bandres
URL:    www.mabandres.com

%}

%% Internal Parameters
a=-a;                                         % to agree with our numerical method
[largo,ancho]=size(z); z=z(:)';      % Input Matrix -> Input Vector
index_paridad=-(z<0)+(z>=0); 
zmax=max(abs(z));
h=0.1;                                       % control number of analytic continuations
nca=ceil(zmax/h);                       % Number of analytic continuations 
[z,index]=sort(abs(z));                 % sort the input vector


%% Analytic Continuation
yi=[1,0];        % Boundary Condition
y=[]; yp=[];    % Initialized variables

for j=1:nca
    xo=(j-1)*h;
    xf=j*h;
    zi=z(xo<=z&z<xf); 
    zi=[zi,xf];
    [yin,ypin]=yc_parabolic(yi,xo,a,zi);
    pint=length(yin)-1;  
    y=[y,yin(1:pint)]; 
    yp=[yp,ypin(1:pint)];       
    yi=[yin(pint+1),ypin(pint+1)]; 
end

if zmax==(nca*h)&&(nca~=0)
    num=sum((z==zmax));
    y=[y,yin(pint+1)*ones(1,num)];
    yp=[yp,ypin(pint+1)*ones(1,num)];
end    
    
if nca~=0
    yv(index)=y;
    ypv(index)=yp;
    ypv=ypv.*index_paridad;
else
    yv=1;
    ypv=0;
end

%% Ouput
y=reshape(yv,[largo,ancho]);     % reshape output to original format
yp=reshape(ypv,[largo,ancho]);
    
return;

function [y,yp]=yO_parabolic(a,z)       % ODD Weber functions (parabolic cylinder functions) 
%{
This function calculates the ODD Weber functions (parabolic cylinder functions) 
which are solutions of the Weber equation (parabolic cylinder equation):

d^2y/dz^2+(z^2/4-a)y=0

where
  -Infinity<z<Infinity 
   a=real parameter

INPUTS:
-Infinity<z<Infinity  REAL independent variable (Vector or Matrix). 
a = real parameter

OUTPUTS:
y =  EVEN Weber function
yo = First derivative of the EVEN Weber function

For more information about Weber functions:

    Miguel A. Bandres and B.M. Rodriguez-Lara, "Nondiffracting accelerating waves: Weber waves and parabolic momentum"
    New Journal of Physics, 15(013054) (2013) ( http://goo.gl/XfaVqq )

    Miguel A. Bandres, J. C. Gutiérrez-Vega, and S. Chávez-Cerda, "Parabolic nondiffracting optical wavefields"
    Optics Letters, 29(1), 44-46 (2004) ( http://goo.gl/KhmqQY )

AUTHOR: Miguel A. Bandres
URL:    www.mabandres.com

%}

%% Internal Parameters
a=-a;                                         % to agree with our numerical method
[largo,ancho]=size(z); z=z(:)';      % Input Matrix -> Input Vector
index_paridad=-(z<0)+(z>=0); 
zmax=max(abs(z));
h=0.1;                                       % control number of analytic continuations
nca=ceil(zmax/h);                       % Number of analytic continuations 
[z,index]=sort(abs(z));                 % sort the input vector


%% Analytic Continuation
yi=[0,1];        % Boundary Condition
y=[]; yp=[];    % Initialized variables

for j=1:nca
    xo=(j-1)*h;
    xf=j*h;
    zi=z(xo<=z&z<xf); 
    zi=[zi,xf];
    [yin,ypin]=yc_parabolic(yi,xo,a,zi);
    pint=length(yin)-1;  
    y=[y,yin(1:pint)]; 
    yp=[yp,ypin(1:pint)];       
    yi=[yin(pint+1),ypin(pint+1)]; 
end

if zmax==(nca*h)&&(nca~=0)
    num=sum((z==zmax));
    y=[y,yin(pint+1)*ones(1,num)];
    yp=[yp,ypin(pint+1)*ones(1,num)];
end    
    
if nca~=0
    yv(index)=y;
    ypv(index)=yp;
    yv=yv.*index_paridad;
else
    yv=0;
    ypv=1;
end

%% Ouput
y=reshape(yv,[largo,ancho]);     % reshape output to original format
yp=reshape(ypv,[largo,ancho]);
    
return;

function [y,yp]=yc_parabolic(yi,xo,a,r) % Subfuction to Calculate Polynomial Expansion
N=100; % Number of polynomial expansion coefficients 

% Calculate coefficients
ac=parabolic_coe(yi,xo,a,N); 
acn=ac./gamma(1:N); 
acp=ac(2:N)./gamma(1:N-1); 

% Evaluate polynomials
y=polyval(acn(N:-1:1),(r-xo));
yp=polyval(acp(N-1:-1:1),(r-xo));
return; 

function [yc]=parabolic_coe(yi,xo,a,N) % Subfuction to Calculate coefficients of Polynomial Expansion 
% Calculate the first N coefficientes of the polynomial expansion

yc=zeros(1,N);
yi=[yi,-((xo^2)/4+a)*yi(1)]; 
yc(1:3)=yi(1:3);
yc(4)=-(4*a+xo^2)/4*yc(2)-xo/2*yc(1);

for n=4:N-1
     yc(n+1)=-(4*a+xo^2)/4*yc(n+1-2)-(n-2)/2*xo*yc(n+1-3)-(n-2)*(n-3)/4*yc(n+1-4);
 end;
 
return;