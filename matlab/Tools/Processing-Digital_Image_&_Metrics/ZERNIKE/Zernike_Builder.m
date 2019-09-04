function [RZS,RMS]=Zernike_Builder(a,A,S,varargin)
% Zernike_Reconstruction is a function that computes 
% a wavefront reconstruction using Zernike's Polynomials
% and by calculating the function expansion coefficients.
%
% Inputs: M, Total number of polynomials to be used
%            WaFr - Square matrix - Wavefront Matrix data
%
% Output, RZS, Reconstructed Zernike Surface,
%            a, expansion coefficients.
%
% Version: 2.0.
%
% Dependencies: zernikes.m - J.J.Cadavid.
%
% Author: Juan Jos?Cadavid Muñoz - Eafit University
% Date - 2014/03/08
%
%% Settings
if nargin == 4
  plotsON = varargin{1};
else
  plotsON = true;
end
%% Initializing
RZS = NaN(S);
M = length(a);
%% Zernike polynomial calculation
[VZk,Pupil] = zernikes(M,A,S);
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
  bar((1:M),a(:));
  axis square;
  title('Zernike Expansion Coefficients','FontSize',20,'FontWeight','bold');
  xlabel('Polynomial number','FontSize',15,'FontWeight','bold');
  ylabel('Amplitude','FontSize',15,'FontWeight','bold');
  set(gca,'fontsize',15);
  grid on
end
end