%% Getfilename
fileName = inputdlg('Input File Name');
SingleFrame = getsnapshot(vid);
imwrite(SingleFrame,fileName{1});