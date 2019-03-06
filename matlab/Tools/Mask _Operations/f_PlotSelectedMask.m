function [mask,wrapMask,wrapMaskFig,maskName] = f_PlotSelectedMask(X,Y, ...
r,phi,gl,phaseValues,tc,s,ph0,p,WsizeRatio,L,f_FR,bcst,period,T0, ...
frkTyp,Aalpha,Angalp,Angbet,z_coeff,z_a,z_pupil,z_disp_wrap, ...
z_plot, normMag,binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx, ...
coordType,abs_ang,MaxMask,plotMask,maskSel)    
% 
% Inputs:
%  Explained inside each function on every case
%  maskSel: selects a specific mask
%
% Outputs:
%  mask: complex structure that has not been truncated and is wrapped on
%        [-pi,pi]. mask = exp(i*UnwrappedMask).
%  wrapMask: truncation, discretization and angle operations on mask. 
%            It is wrapped on [-pi,pi].
%  wrapMaskFig: figure handler if needed outside the function

%  maskName: name of the selected mask

%% Mask selection
switch maskSel 
 case 0 % Spiral phase mask or mapa de fase espiral or máscara espiral
        % de fase or máscara helicoidal de fase 
  [mask,wrapMask,wrapMaskFig] = f_SpiralMask(r,phi,gl,phaseValues,tc,s, ...
  ph0,normMag,binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx, ...
  coordType,abs_ang,MaxMask,plotMask);
  maskName = 'Spiral';
   
 case 1 % Laguerre-Gauss (LG) beams
  [mask,wrapMask,wrapMaskFig] = f_LGMask(r,phi,gl,phaseValues,tc,s,ph0, ...
  p,WsizeRatio,normMag,binMask,binv,MaskPupil,rSize,monitorSize, ...
  scrnIdx,coordType,abs_ang,MaxMask,plotMask);
  maskName = 'LG';
 
 case 2 % Vortex Producing Lens (VPL) = Helicoidal + Fresnel lens
  [mask,wrapMask,wrapMaskFig] = f_VPLMask(r,phi,gl,phaseValues,tc,s, ...
  ph0,L,f_FR,normMag,binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx, ...
  coordType,abs_ang,MaxMask,plotMask);
  maskName = 'VPL';
  
 case 3 % Elliptic Gaussian Vortex (EGV) mask or Vórtice elíptico-gaussiano
  [mask,wrapMask,wrapMaskFig] = f_EGVMask(X,Y,r,gl,phaseValues,tc,s, ...
  ph0,bcst,normMag,binMask,binv,MaskPupil,rSize,monitorSize,scrnIdx, ...
  coordType,abs_ang,MaxMask,plotMask);
  maskName = 'EGV';
  
 case 4 % Threefold dislocation hologram or double pitch fork hologram or
        % fork phase mask or holograma en forma de tenedor or fork
        % grating plate or rejilla de Ronchi con una dislocación 
        % (2013_Vortex_Generations_Mach-Zehnder_Interferometer)
  [mask,wrapMask,wrapMaskFig] = f_ForkMask(X,Y,r,phi,gl,phaseValues,tc, ...
  s,ph0,period,T0,frkTyp,Aalpha,Angalp,Angbet,normMag,binMask,binv, ...
  MaskPupil,rSize,monitorSize,scrnIdx,coordType,abs_ang,MaxMask,plotMask);
  maskName = 'Fork';
  
%%%%%%%%%%%%%%%%%%%% NOT USED BUT FOR ACADEMIC PURPOSES %%%%%%%%%%%%%%%%%%%
 case 5 % Zernike (aberrations)
  [mask,wrapMask,wrapMaskFig] = f_ZernikeMask(X,Y,r,z_coeff,z_a,L,gl, ...
  phaseValues,z_pupil,z_disp_wrap,z_plot,normMag,binMask,binv,rSize, ...
  monitorSize,scrnIdx,coordType,MaxMask,plotMask);
  maskName = 'Zernike';
                  
 case 6 % Laguerre-Gauss (LG) + Zernike
  [mask,wrapMask,wrapMaskFig] = f_LGZernikeMask(X,Y,r,phi,gl, ...
  phaseValues,tc,s,ph0,p,WsizeRatio,z_coeff,z_a,L,z_pupil, ...
  z_disp_wrap,z_plot,normMag,binMask,binv,rSize,monitorSize,scrnIdx, ...
  coordType,abs_ang,MaxMask,plotMask);
  maskName = 'LG_Zernike';       
  
 case 7 % Hermite-Gauss (HG) beams
  % maskName = 'HG'; 
  error('Mask not available yet');

 case 8 % Mutliple vortices
  % Maybe use mirror padarray!
  % There would be superposition or one may need two SLM's
  % maskName = 'MultipleVortx'; 
  error('Mask not available yet');
  
 case 9 % Sum of spiral phase masks
  % Check: 1_2013_IMP_Storing_High-Dimensional_Quantum_States_in_a_Cold ...
  %        Atomic Ensemble_Ding
  % maskName = 'SpiralPhaseSummed'; 
  error('Mask not available yet');
  
 case 10 % Gerchberg-Saxton
  % Function itself is done, add and edit constraints 
  % Gerchberg_Saxton;   
  % maskName = 'GS'; 
  error('Mask not available yet');    
  
 otherwise % Void; Ideal response; no aberrations
  mask = 1; % Unitary
  maskName = 'UnitaryFilter';
end
end