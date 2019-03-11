function [RZS,RMS] = f_ZernikeBuilder(X,Y,a,apperture,varargin)
% Zernike_Reconstruction is a function that computes 
% a wavefront reconstruction using Zernike's Polynomials
% and by calculating the function expansion coefficients.
%
% Inputs: 
%  X,Y: A grid of the spatial vector: 2D Cartesian coordiantes. They
%  may already have the shiftCart
%  a: Zernike weight vector that contains the aberrations contributions
%  apperture: float that defines pupil relative size (w.r.t. sSize), like a
%     percentage
%  varargin: (0) no plots; (1) plots
% 
% Outputs:
%         RZS - Reconstructed Zernike Surface,
%         RMS - Root Mean Square
%
% Version: 2.0.
%
% Dependencies: f_zernikes.m - J.J.Cadavid.
%
% Author: Juan Jose Cadavid Munoz - Eafit University
% Date - 2014/03/08
%
% Commented by Samuel Plazas Escudero on 2018/04/03 and variables
% X,Y were added on 2019/02/27

%% Settings
if nargin == 5
  plotsON = varargin{1};
else
  plotsON = true;
end

%% Initializing
% RZS = NaN(spaceSize); % Initialization: SxS matrix
sizX = size(X); % Equals Y's.
spatialSize = min(sizX);  % Square matrix of the total size of the figure
                          % Minimum since that's the size of a square that
                          % can fit inside the screen's rectangle
RZS = zeros(spatialSize); % Initialization: SxS matrix
p = length(a); % MxM matrix will be created
               % Correction in order to get correctly the polynomials

%% Zernike polynomial calculation
[VZk,Pupil] = f_ZernikePolynomials(X,Y,p,apperture); 
% Creates Z. polynomials

%% Wavefront Reconstruction
RZSV = sum(bsxfun(@times, a.', VZk), 2); 
RZS(~isnan(Pupil)) = RZSV;

%% RMS Calcuation
RMS = sqrt(sum(a .^ 2));

%% Plotting
 if plotsON
  figure;
  imagesc(RZS),colorbar;
  axis square;
  title('Zernike Reconstructed Surface','FontSize',20,'FontWeight','bold');
  xlabel('X Deformations - Pixels','FontSize',15,'FontWeight','bold');
  ylabel('Y Deformations - Pixels','FontSize',15,'FontWeight','bold');
  set(gca,'fontsize',15);

  figure;
  bar((0:p-1),a(:)); % -1 corrects a deviation of one of the coefficients
  axis square;
  title('Zernike Expansion Coefficients','FontSize',20,'FontWeight','bold');
  xlabel('Polynomial number','FontSize',15,'FontWeight','bold');
  ylabel('Amplitude','FontSize',15,'FontWeight','bold');
  set(gca,'fontsize',15); xlim([0 p-1])
  grid on
 end
end