% clear all; close all; clc
%%%%%%%%%%%%%%%%%%%%% This code makes uniform swaths along a line created from two points that the user selects on a map

numberofswaths = 4; %set the number of swaths you want to make

%% Initiate some stuff
COLORMAPtopo =colormap_cpt('Topo'); % colormap, 'Topo' is a cpt file in my directory

%% Open DEM and pick points to plot swaths
imageschs(DEM, [],'colormap',COLORMAPtopo,...
    'colorbar', true,'caxis',[0 1400])
hold on;

[coordsX, coordsY] = ginput(2);
upleft_y = 2.196e6;
%%
% You can manually enter coordinates here... if you want...
% coordsX(1) = 283272;
% coordsY(1) = 2174845;
% coordsX(2) = 344564
% coordsY(2) = 2147163;
uplefty = 2.196e6; 
x = coordsX(1):1000:coordsX(2);
x1 = x - coordsY(1);
slp = (coordsY(2)-coordsY(1))/(coordsX(2)-coordsX(1));
y = (slp*x1)+coordsY(1);
% dist = sqrt((coordsX(2)-coordsX(1))^2 + (coordsY(2) - coordsY(1))^2);
xspace = (coordsX(2) - coordsX(1))/numberofswaths;
xdiv = coordsX(1):xspace:coordsX(2);
xdiv1 = xdiv - coordsX(1);
ydiv = (slp*xdiv1)+coordsY(1);
plusy = uplefty-coordsY(1);
newx = [];
newy = [];
for i = 1:length(xdiv)
    newy(i) = ydiv(i) + plusy;
    newx(i) = (plusy*(-slp))+xdiv(i);
end
imagesc(DEM)
hold on; plot(x,y, 'k', 'LineWidth', 5)
plot(xdiv, ydiv, 'or')
plot(newx, newy, '*')

%%
SWATHS = struct;
for i = 1:length(newx)
    SWATHS.ind(i) = SWATHobj(DEM, [xdiv(i),newx(i)], [ydiv(i),newy(i)])
end
%%
close all
SW = struct;
swidth = 1e4;
vexag = 5;
for i = 1:length(SWATHS.ind)
    [SWATH,SwathMats,XYpoints,Bends]=MakeTopoSwath...
        (DEM,SWATHS.ind(i).xy0,swidth,...
        'plot_figure',true,'vex',vexag, 'plot_as_points', true)
    xlabel = ('Distance (m)')
    ylabel = ('Elevation (m)')
    
    SW.swath(i) = SWATH;
end
%%
figure(6)
COLORMAPtopo =colormap_cpt('Topo');
imageschs(DEM, [],'colormap',COLORMAPtopo,...
    'colorbar', true,'caxis',[0 1400])
hold on; 
for i = 1:length(SWATHS.ind)
    SW.swath(i).distx = SW.swath(i).distx./1000;
    A = SW.swath(i)
    h = plot(A,'outline',true)
    for ind = 1:3
        h(ind).Color = 'k';
        h(ind).LineWidth = 1;
    end  
end
