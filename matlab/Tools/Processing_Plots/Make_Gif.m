function Make_Gif(outfile,t,delay)
%% gif utilities
set(gcf,'color','w'); % set figure background to white
drawnow;
frame = getframe(1);
im = frame2im(frame);
[imind,cm] = rgb2ind(im,256);
% outfile = 'sinewave.gif';
%% On the first loop, create the file. In subsequent loops, append.
if t==1
    imwrite(imind,cm,outfile,'gif','DelayTime',delay,'loopcount',inf);
else
    imwrite(imind,cm,outfile,'gif','DelayTime',delay,'writemode','append');
end
%EOF
end