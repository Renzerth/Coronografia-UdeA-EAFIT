function [RZS,a,RMS,Error] = f_Zernike_Reconstruction(M,WaFr,A)
% Zernike_Reconstruction is a function that computes 
% a wavefront reconstruction using Zernike's Polynomials
% and by calculating the function expansion coefficients.
%
% Inputs: M, Total number of polynomials to be used
%         WaFr - Square matrix - Wavefront Matrix data
%         A - float - Defines pupil relative size (w.r.t. S)
%
% Output, RZS, Reconstructed Zernike Surface,
%         a, expansion coefficients.
%
% Version: 2.0.
%
% Dependencies: zernikes.m - J.J.Cadavid.
%
% Author: Juan José Cadavid Muñoz - Eafit University
% Date - 2014/03/08
% Commented by Samuel Plazas Escudero on 2018/04/03

%% Initializing
S = size(WaFr,1);
RZS = NaN(S);
% WaFrP = zeros(S);

%% Zernike polynomial calculation
[VZk,Pupil] = f_zernikes(M,A,S);

%% Definition of Wavefront in the Pupil
WaFr = WaFr + Pupil;

%% Matrix-Vector convertion & NaN deletion
WaFrV = WaFr(~isnan(WaFr));

%% Coefficient decomposition - SVD
a = VZk\WaFrV;

%% Wavefront Reconstruction
RZSV = sum(bsxfun(@times, a.', VZk), 2);
RZS(~isnan(Pupil)) = RZSV;
% for i=4:M % Reconstruction with piston or tilts
%     RZS(~isnan(Pupil)) = VZk(:,i);
%     WaFrP=WaFrP+a(i)*RZS;
% end

%% RMS Calcuation
RMS = sqrt(sum(a .^ 2));

%% Normalization of Deformations
NRZS = RZS-min(min(RZS));
NRZS = NRZS/max(max(NRZS));
NWaFr = WaFr-min(min(WaFr));
NWaFr = NWaFr/max(max(NWaFr));

%% Absolute error of reconstruction
Error = abs(NWaFr-NRZS);

%% Ploting
figure
imagesc(RZS),colorbar;
axis square;
title('Zernike Reconstructed Surface','FontSize',20,'FontWeight','bold');
xlabel('X Deformations - Pixels','FontSize',15,'FontWeight','bold');
ylabel('Y Deformations - Pixels','FontSize',15,'FontWeight','bold');
set(gca,'fontsize',15);

figure
bar((1:M),a(:));
axis square;
title('Zernike Expansion Coefficients','FontSize',20,'FontWeight','bold');
xlabel('Polynomial number','FontSize',15,'FontWeight','bold');
ylabel('Amplitude','FontSize',15,'FontWeight','bold');
set(gca,'fontsize',15);
grid on

figure
imagesc(Error),colorbar;
axis square;
title('Absolute Reconstruction Error','FontSize',20,'FontWeight','bold');
xlabel('X Deformations - Pixels','FontSize',15,'FontWeight','bold');
ylabel('Y Deformations - Pixels','FontSize',15,'FontWeight','bold');
set(gca,'fontsize',15);
end