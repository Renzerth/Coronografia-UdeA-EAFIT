function f_ComputeMaskGradient(x,y,mask,gradMask)
% computes the gradient of a real-valued matrix (mask)
% Plots the net gradient and the x, y profiles
% Inputs: x,y coordinates vectors, Real-valued matrix and boolean for
% activating the function (gradMask)
if gradMask == 1
 [xg,yg] = gradient(mask);
 figure; contour(x,y,mask); hold on;
 quiver(x,y,xg,yg); title('Gradient of the mask'); hold off; colorbar;
 figure; imagesc(x,y,xg); 
 colormap(hot); title('X-profile gradient of the mask');
 figure; imagesc(x,y,yg); 
 colormap(hot); title('Y-profile gradient of the mask');
end
end