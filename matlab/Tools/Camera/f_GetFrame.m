% like preview
function [Frame] = f_GetFrame(vid)
preview(vid);
start(vid);
stoppreview(vid);
closepreview(vid);
Frame = getdata(vid);
end