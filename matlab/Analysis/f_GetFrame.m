function [Frame] = f_GetFrame(vid)
preview(vid);
start(vid);
stoppreview(vid);
Frame = getdata(vid);
end