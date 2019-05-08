function [rearranged] = flipLinearIndexCell(cellData,width)
cellLength = length(cellData);
rearranged = cell(1,cellLength);

for widthIndex = 1:width
    for dataIndex = 1:width
        flipDataIndex = widthIndex + width*fix(dataIndex/width);
        arrangeIndex = dataIndex + width*(widthIndex-1);
        rearranged{arrangeIndex} = cellData{flipDataIndex};
    end
    dataGroup = 6*fix(widthIndex/6);
end