clear;close all;clc;

proj = projcrs(32610);

M = dlmread('../meshes/Turning_River_Dam_Burn_1150.bcode');
M = M(2:end,:);
load('../Turning_DEM.mat');

% Read dam shapefile
dam1 = shaperead('../dams/dam1.shp');
dam1.X = dam1.X + 0.00025;
dam1.Y = dam1.Y - 0.001;
dam2 = shaperead('../dams/dam2.shp');
dam2.X = dam2.X + 0.00036;
dam2.Y = dam2.Y - 0.00071;
% Project dam coordinates
[dam1.X,dam1.Y] = projfwd(proj,dam1.Y,dam1.X);
[dam2.X,dam2.Y] = projfwd(proj,dam2.Y,dam2.X);
ind = find(dam2.Y < 3.627*1e6);
dam2.X(ind) = []; dam2.Y(ind) = [];

d1 = polybuffer([dam1.X' dam1.Y'],'lines',50);
d2 = polybuffer([dam2.X' dam2.Y'],'lines',50);

figure;
plot(M(:,1),M(:,2),'kx'); hold on;
%plot(dam1.X,dam1.Y,'r-','LineWidth',2);  plot(dam2.X,dam2.Y,'g-','LineWidth',2); 
plot(d1.Vertices(:,1),d1.Vertices(:,2),'r-','LineWidth',2);
plot(d2.Vertices(:,1),d2.Vertices(:,2),'g-','LineWidth',2);

in1 = inpoly2([xdem(:) ydem(:)],d1.Vertices);
in2 = inpoly2([xdem(:) ydem(:)],d2.Vertices);
in = in1 | in2;


figure;
scatter(xdem(in),ydem(in),18,dem(in),'filled'); colorbar; hold on;

in1 = inpoly2([M(:,1) M(:,2)],d1.Vertices);
in2 = inpoly2([M(:,1) M(:,2)],d2.Vertices);
in = in1 | in2;
in = find(in == 1);

figure;
scatter(M(in,1),M(in,2),18,M(in,3),'filled'); colorbar; hold on; clim([24 36]);

for i = 1 : length(in)
    pt = polybuffer([M(in(i),1) M(in(i),2)],'points',50);
    %plot(pt.Vertices(:,1),pt.Vertices(:,2),'r-','LineWidth',1); hold on;
    indam = inpoly2([xdem(:) ydem(:)],pt.Vertices);
    M(in(i),3) = max(dem(indam))*1.5;
end


figure;
scatter(M(in,1),M(in,2),18,M(in,3),'filled'); colorbar; hold on; clim([24 36]);

M(:,3) = round(M(:,3),1);

fout = '../meshes/Turning_River_Dam_Burn_1150';
num_of_nodes = size(M,1);
fprintf('\n Writing to .bcode file...');
fileID = fopen(strcat(fout,'.bcode'),'w');
fprintf(fileID, '%d \n', num_of_nodes);
fprintf(fileID, '%.2f \t %.2f \t %.2f \t %d \n', [M(:,1) M(:,2) M(:,3) M(:,4)]');
fclose(fileID); 
fprintf('Done!\n');
    
