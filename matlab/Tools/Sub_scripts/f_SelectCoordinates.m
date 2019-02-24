function [X,Y,r,phi] = f_SelectCoordinates(Xslm,Yslm,rSLM,phiSLM,Xpc,...
                                           Ypc,rPC,phiPC,plotMask)
%% Coordinates selection
if plotMask == 2 % SLM
    r = rSLM; phi = phiSLM;
    X = Xslm; Y = Yslm;
else % plotMask == 0 or 1 or 3 % PC
    r = rPC; phi = phiPC; 
    X = Xpc; Y = Ypc;    
end

end

