function [rearranged] = flipLinearIndexCell(cellData, rowWidth, columnWidth)
rearranged = cell(1,length(cellData));

for rowIndex =1:rowWidth
    for columnIndex = 1:columnWidth
        dataIndex = rowIndex + (columnIndex-1)*rowWidth;
        arrangingtIndex = columnIndex + (rowIndex-1)*columnWidth;
        rearranged{arrangingtIndex} = cellData{dataIndex};
    end
end
end