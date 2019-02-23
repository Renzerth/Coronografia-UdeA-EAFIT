switch maskSel 
    
 case 0 % Spiral phase mask or mapa de fase espiral or m�scara espiral
        % de fase or m�scara helicoidal de fase 
  mask = f_SpiralMask(r,phi,gl,glphi,mingl,maxgl,levShft,tc,s,ph0, ...
                       binMask,showM);
  maskName = 'Spiral';
   
 case 1 % Laguerre-Gauss (LG) beams
  mask = f_LGMask(r,phi,gl,glphi,mingl,maxgl,levShft,tc,s,ph0,p,W, ...
                   binv,norm,abs_ang,binMask,monitorSize,scrnIdx,showM);
  maskName = 'LG';
  
 case 2 % Vortex Producing Lens (VPL) = Helicoidal + Fresnel lens
  mask = f_VPLMask(r,phi,gl,glphi,mingl,maxgl,levShft,tc,s,ph0,L, ...
                    f_FR,binMask,showM);
  maskName = 'VPL';
  
 case 3 % Elliptic Gaussian Vortex (EGV) mask or V�rtice el�ptico-gaussiano
  mask = f_EGVMask(X,Y,r,gl,glphi,mingl,maxgl,levShft,tc,s,ph0, ...
                    bcst,binMask,showM);
  maskName = 'EGV';
  
 case 4 % Threefold dislocation hologram or double pitch fork hologram or
        % fork phase mask or holograma en forma de tenedor or fork
        % grating plate or rejilla de Ronchi con una dislocaci�n 
        % (2013_Vortex_Generations_Mach-Zehnder_Interferometer)
  mask = f_ForkMask(X,Y,r,phi,gl,glphi,mingl,maxgl,levShft,tc,s, ...
                ph0,L,period,T0,frkTyp,Aalpha,Angalp,Angbet,binMask,showM);
  maskName = 'Fork';
  
%%%%%%%%%%%%%%%%%%%% NOT USED BUT FOR ACADEMIC PURPOSES %%%%%%%%%%%%%%%%%%%
 case 5 % Zernike (aberrations)
  sSize = min(min(size(X)),min(size(Y)));
  mask = f_ZernikeMask(x,y,r,z_coeff,a,frac,L,gl,glphi,mingl,maxgl, ...
                      levShft,pupil,sSize,disp_wrap,plot_z, ...
                      binMask,monitorSize,scrnIdx,showM);
  
 case 6 % Laguerre-Gauss (LG) + Zernike
  sSize = min(max(X),max(Y));
  mask = f_LGZernikeMask(x,y,r,phi,gl,glphi,mingl,maxgl,levShft,tc,s, ...
                           ph0,p,W,binv,norm,abs_ang,z_coeff,a,frac,L, ...
                           pupil,sSize,disp_wrap,plot_z,binMask, ...
                           monitorSize,showM);
                       
 case 7 % Hermite-Gauss (HG) beams
     
 case 8 % Mutliple vortices
  % Maybe use mirror padarray!
  % There would be superposition or one may need two SLM's
 
 case 9 % Sum of spiral phase masks
  % Check: Storing_High-Dimensional_Quantum_States_in_a_Cold_
  
 case 10 % Gerchberg-Saxton
  % Function itself is done, add and edit constraints 
  % Gerchberg_Saxton;   
                       
 otherwise % Void; Ideal response; no aberrations
  mask = 1; % Unitary
  
end