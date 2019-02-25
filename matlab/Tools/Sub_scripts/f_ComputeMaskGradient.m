function f_ComputeMaskGradient(x,y,mask,gradMask)
if gradMask == 1
 [xg,yg] = gradient(angle(mask));
 figure; contour(x,y,angle(mask)); hold on;
 quiver(x,y,xg,yg); title('Gradient of the mask'); hold off; colorbar;
 figure; imagesc(x,y,xg); 
 colormap(hot); title('X-profile gradient of the mask');
 figure; imagesc(x,y,yg); 
 colormap(hot); title('Y-profile gradient of the mask');
end
end