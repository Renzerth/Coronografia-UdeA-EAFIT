function regionCentroid = f_approximateSpotSize(cdata)
%% Estimación del tamaño en pixeles del disco de Airy
binaryData = im2bw(cdata,0.01);
regionInfo = regionprops(binaryData,'Centroid','area','MajorAxisLength','MinorAxisLength');
[~, sortIndexes] = sort(cat(1,regionInfo.Area), 'descend');
mainIndex = sortIndexes(1);
blobRadius = mean([regionInfo(mainIndex).MajorAxisLength, regionInfo(mainIndex).MinorAxisLength],2)/2;
regionCentroid = cat(1,regionInfo(mainIndex).Centroid);
radiusArea = sqrt(max(regionInfo(mainIndex).Area)/pi); 
aproxRadius = ceil(mean([blobRadius,radiusArea])); % Once a pixel is occupied, it is assumed to fill the whole pixel
%% Proporción disco de Airy - tamaño de imagen
[p,q] = size(cdata);
halfP = p/2;
halfQ = q/2;
xLengthRatio = halfP/aproxRadius;
yLengthRatio = halfQ/aproxRadius;
%% Coordenadas escaladas
x = -xLengthRatio:xLengthRatio-1;
y = -yLengthRatio:yLengthRatio-1;
%% plot
imagesc(x,y,cdata);
end