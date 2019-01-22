%% Example of the IG.Ince_Gaussian function
%{
This script plots the first 25 even/odd Ince Gaussian Beams 
using the IG.Ince_Gaussian function:

[IGB,X,Y]=IG.Ince_Gaussian(L,N,parity,p,m,e,w0,k,z)

%}

close all;

%% Ellipticity 
ec=2;  % Change the ellipticity to see how the Ince Gaussian Modes behave, for example:
       % ec=0.0001; if ec->0   then Ince Gaussian Beams -> Laguerre Gaussian Beams
       % ec=1000;   if ec->Inf then Ince Gaussian Beams -> Hermite Gaussian Beams

%% First 25 EVEN Ince Gaussian Beams 
figure; axis off;
title('EVEN Ince Gaussian Modes','FontSize',12);
count=1; 
for p=0:7
    for m=(0:2:p)+mod(p,2)
        [IGB,X,Y]=f_InceGaussian(15e-3,201,0,p,m,ec,3e-3,(2*pi/632.8e-9),0); 
        hax=axes();
        surf(X,Y,abs(IGB)); shading interp; lighting phong; view(2); axis equal; axis tight;  axis off;
        rowIdx=fix((count-1)/5); colIdx=mod(count-1,5); newPos=[.02+.16*colIdx,0.7-.2*rowIdx,.25,.25]; set(gca,'outer',newPos),
        count=count+1;
    end
end

%% First 25 ODD Ince Gaussian Beams 
figure; axis off;
title('ODD Ince Gaussian Modes','FontSize',12);
count=1;
for p=1:8
    for m=(2:2:(p+mod(p,2)*2))-mod(p,2)
        [IGB,X,Y]=f_InceGaussian(15e-3,201,1,p,m,ec,3e-3,(2*pi/632.8e-9),0); 
        hax=axes();
        surf(X,Y,abs(IGB)); shading interp; lighting phong; view(2); axis equal; axis tight;  axis off;
        rowIdx=fix((count-1)/5); colIdx=mod(count-1,5); newPos=[.02+.16*colIdx,0.7-.2*rowIdx,.25,.25]; set(gca,'outer',newPos),
        count=count+1;
    end
end