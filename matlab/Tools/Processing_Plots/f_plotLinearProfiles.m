function f_plotLinearProfiles(dataArray,x,y,cartcoordX,cartcoordY,tit,xlab,ylab,plotData,profileData, plotH,plotV,tol,refPoints,fontSize,lineWidth)
%% Plot lines reference
Hx = refPoints(1,:); % Horizontal x components: (xi,xf)
Vy = refPoints(2,:);  % Vertical y components: (yi,yf)
Hy = refPoints(3,:); % Horizontal y components: (yi,yf)
Vx = refPoints(4,:); % Vertical x components: (xi,xf)

%% Profiles of the 2D array drawn on the dataArray
if plotData
    figure; imagesc(x,y,dataArray); title(tit); % colormap(hot)
    xlabel(xlab);
    ylabel(ylab);
    hold on
    % Horizontal:
    line(x(Hx),y(Hy),'LineWidth',lineWidth,'Color','blue','LineStyle','--');
    % Vertical:
    line(x(Vx),y(Vy),'LineWidth',lineWidth,'Color','red','LineStyle','--');
    hold off
end

%% Plot of each profiles
if plotH
    %% Horizontal profile
    figure();
    plot(cartcoordX(1+tol:end-tol),profileData(1+tol:end-tol));
    title(strcat('Horizontal profile of the',{' '},tit),'FontSize',fontSize,'FontWeight','bold');
end

if plotV
    %% Vertical profile
    figure();
    plot(cartcoordY(1+tol:end-tol),profileData(1+tol:end-tol)');
    title(strcat('Vertical profile of the',{' '},tit),'FontSize',fontSize,'FontWeight','bold');
end

end