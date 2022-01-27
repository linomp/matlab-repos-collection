clear all; close all; clc
%[X,Y,Z] = cylinder(10,30);
peaks
%surfc(X,Y,Z,'linestyle','none')
alpha(0.5)
colorbar
colormap(jet)
set(gcf,'color','k')
set(gca,'color','k','YColor','w','XColor','w','ZColor','w') 
C = colorbar;
set(C,'Color','w')
ylabel('y','Color','w')
xlabel('x','Color','w')
zlabel('z','Color','w')