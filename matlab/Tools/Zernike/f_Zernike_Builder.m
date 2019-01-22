function [RZS,RMS] = f_Zernike_Builder(a,A,sSize,varargin)
% Zernike_Reconstruction is a function that computes 
% a wavefront reconstruction using Zernike's Polynomials
% and by calculating the function expansion coefficients.
%
% Inputs: 
%        a - Zernike weight vect: contains the aberrations contributions
%        A - float - Defines pupil relative size (w.r.t. sSize), like a
%                    percentage
%        sSize - Square matrix: total size of the figure. Space Size
%        varargin - 0: no plots; 1: plots
%             
% Outputs:
%         RZS - Reconstructed Zernike Surface,
%         RMS - Root Mean Square
%
% Version: 2.0.
%
% Dependencies: f_zernikes.m - J.J.Cadavid.
%
% Author: Juan Jos� Cadavid Mu�oz - Eafit University
% Date - 2014/03/08
% Commented by Samuel Plazas Escudero on 2018/04/03

%% Settings
if nargin == 4
  plotsON = varargin{1};
else
  plotsON = true;
end

%% Initializing
% RZS = NaN(spaceSize); % Initialization: SxS matrix
RZS = zeros(sSize); % Initialization: SxS matrix
M = length(a); % MxM matrix will be created
                 % Correction in order to get correctly the polynomials

%% Zernike polynomial calculation
[VZk,Pupil] = f_zernikes(M,A,sSize); % Creates Z. polynomials

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
  bar((0:M-1),a(:)); % -1 corrects a deviation of one of the coefficients
  axis square;
  title('Zernike Expansion Coefficients','FontSize',20,'FontWeight','bold');
  xlabel('Polynomial number','FontSize',15,'FontWeight','bold');
  ylabel('Amplitude','FontSize',15,'FontWeight','bold');
  set(gca,'fontsize',15); xlim([0 M-1])
  grid on
 end
end