function f_plotLinearProfiles(dataArray,x,y,tit,xlab,ylab,plotData,plotH,plotV,tol)
%% Profiles of the 2D array drawn on the dataArray
if plotData
    figure; imagesc(x,y,dataArray); title(tit); % colormap(hot)
    xlabel(xlab);
    ylabel(ylab);
    hold on
    % Horizontal:
    line(x(Hx),y(Hy),'LineWidth',3,'Color','blue','LineStyle','--');
    % Vertical:
    line(x(Vx),y(Vy),'LineWidth',3,'Color','red','LineStyle','--');
    hold off
end

%% Plot of each profiles
if plotH
    %% Horizontal profile
    figure();
    plot(x(1+tol:end-tol),Hprof(1+tol:end-tol));
    title(strcat('Horizontal profile of the',{' '},tit));
end

if plotV
    %% Vertical profile
    figure();
    plot(y(1+tol:end-tol),Vprof(1+tol:end-tol));
    title(strcat('Vertical profile of the',{' '},tit));
end

end