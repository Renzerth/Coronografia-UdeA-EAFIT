% Illustrates the use of LaguerreGaussianE.m by plotting the 
% first 16 Laguerre-Gaussian modes (p,m)=([0:3],[0:3]).
%% Parameters
var = 2; % 1: Magnitude; 2: Phase
modified = 0; % 1: adds the same polinomial but with 1 index higher and
              % rotated 90 degrees with a i
p = 3;
m = 3;

%% Sampling and creation of vectors
pnumbers = 0:p-1;
mnumbers = 0:m-1;
pmax=max(pnumbers);
mmax=max(mnumbers);

screensize = 0.08; % 0.08
npts = 75; % Number of samples % 75
domain = zeros(npts,npts,2); % Variable initialization
x=-screensize:2*screensize/(npts-1):screensize;
y=x;
w = 0.022; % 0.022
R=-1e3; % -1e3
deltax=0*w;
deltay=0*w;
wfactor = 1; % 1
lambda=1.064e-6;
ploton=[1,1];
n=40; % 40


[xmesh,ymesh]=meshgrid(x,y);
domain(:,:,1)=xmesh; domain(:,:,2)=ymesh;
for p = pnumbers
for m = mnumbers

subplot(pmax+1,mmax+1,p*(mmax+1)+m+1)
hg = f_HermiteGaussianE([p,m,f_axial_q_(w*wfactor,R,lambda),lambda],xmesh+deltax,ymesh+deltay);

if modified == 1
    hg = hg + 1i.*f_HermiteGaussianE([p+1,m+1,f_axial_q_(w*wfactor,R,lambda),lambda],xmesh+deltax,ymesh+deltay);
end

if var == 1
    h=pcolor(x,y,abs(hg).^2);
else % var == 1
    h=pcolor(x,y,angle(hg));
end
colormap('bone')
set(h,'EdgeColor','none');
set(h,'FaceColor','interp');
set(gca,'Visible','off');
set(gcf,'Color','black');
axis square
shg;

end
end