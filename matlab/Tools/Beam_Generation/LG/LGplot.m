%% Plots Laguerre-Gauss or Optical Vortex beams
%  p: degree of polynomial. Is known as the number of radial nodes in the
%     intensity distribution. There are p+1 rings or airy disks present. 
%     Denotes the rings transitions
%  m: order of the derivative that defines it; positive for the radial part
%     but the sign can vary for the azimuthal: direction of the spiral 
%     mask. Correspons to the topological charge of the spiral phase mask
%     or the phase discontinuity. It is known as the azimuthal index giving
%     an OAM of l*planck_reduced per photon. Vorticity increases with it.

% clc; clear; close all; % Initialization

%% Space Set-Up
S = 512; % Number of samples
x = (-1:(2/(S-1)):1); % Makes an array with the size of the WaFr matrix
[X,Y] = meshgrid(x,x);
[phi,r] = cart2pol(X,Y);

%% Parameters
W = 150; % Width of the modes: related with the radius of the phase and on
         % the disks on the magnitude
W = W/S; % Normalization with number of samples
ph0 = 0; % Initial phase of the angle (or spiral phase mask)
s = +1; % Sign of the spiral phase mask (+1 or -1)
modes = 2; % Modes of polynomials
var = 2; % 1: Magnitude; 2: Phase

%% Calculation
figure(var);
for p = 0 : modes % y-axis variation
    for m = 0 : modes % x-axis variation;
        subplot(modes+1,modes+1,p*(modes+1)+m+1)
        if var == 1
            h = pcolor(X,Y,angle(f_LaguerreGauss(r,phi,m,s,ph0,p,W)));
            % Wrapped on [-pi,pi]; modulo(2pi)
        else % var == 1
            h = pcolor(X,Y,abs(f_LaguerreGauss(r,phi,m,s,ph0,p,W)));
            % Only positive values since it is a intensity measurement
        end
        colormap('bone'), set(h,'EdgeColor','none');  
        set(h,'FaceColor','interp'), set(gca,'Visible','off');
        set(gcf,'Color','black'), axis square, hold off; shg;
        % shg makes the current figure visible and raises it above all 
        % other figures on the screen.
    end
end