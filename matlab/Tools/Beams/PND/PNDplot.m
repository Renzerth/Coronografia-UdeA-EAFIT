%% Plot PND
[PB,X,Y] = f_Parabolic_Nondiffracting_Beam(6*10^-3,501,0,10^4,3); figure;
surf(X,Y,abs(PB)); shading interp; lighting phong; view(2); colormap(hot);
axis equal tight; xlabel('x','FontName','Helvetica','Fontsize',15); 
ylabel('y','FontName','Helvetica','Fontsize',15); 
ht=title('Parabolic Nondiffracting Beam','Fontsize',12);
set(ht,'FontName','Helvetica','FontSize',12,'FontWeight','bold'); 
set(gca,'FontName','Helvetica','FontSize',10,'FontWeight','bold'); 
