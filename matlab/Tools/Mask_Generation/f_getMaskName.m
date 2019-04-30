function [maskName] = f_getMaskName(maskSel)
%% Mask selection
switch maskSel
    case 0 % Spiral phase mask or mapa de fase espiral or mascara espiral
        maskName = 'Spiral';
        
    case 1 % Laguerre-Gauss (LG) beams
        maskName = 'LG';
        
    case 2 % Vortex Producing Lens (VPL) = Helicoidal + Fresnel lens
        maskName = 'VPL';
        
    case 3 % Elliptic Gaussian Vortex (EGV) mask or Vortice eliptico-gaussiano
        maskName = 'EGV';
        
    case 4 % Threefold dislocation hologram or double pitch fork hologram
        maskName = 'Fork';
        
    case 5 % Zernike (aberrations)
        maskName = 'Zernike';
        
    case 6 % Laguerre-Gauss (LG) + Zernike
        maskName = 'LG_Zernike';
        
    case 7 % Hermite-Gauss (HG) beams
        error('Mask not available yet');
        
    case 8 % Mutliple vortices
        error('Mask not available yet');
        
    case 9 % Sum of spiral phase masks
        error('Mask not available yet');
        
    case 10 % Gerchberg-Saxton
        error('Mask not available yet');
        
    otherwise % Void; Ideal response; no aberrations
        maskName = 'UnitaryFilter';
end
end