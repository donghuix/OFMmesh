clear;close all;clc;

addpath('/Users/xudo627/Developments/mylib/m/');


vals = [11;12;21;22;23;24;31;41;42;43;51;52;71;72;81;82;90;95];
class =  {'Open Water', 'Perennial Ice and Snow','Developed, Open Space', ...
          'Developed, Low Intensity', 'Developed, Medium Intensity', ...
          'Developed, High Intensity', 'Barren Land (Rock/Sand/Clay)', 'Deciduous Forest', ...
          'Evergreen Forest', 'Mixed Forest', 'Dwarf Scrub', 'Shrub/Scrub', 'Grassland/Herbaceous', ...
          'Sedge/Herbaceous', 'Pasture/Hay', 'Cultivated Crops', 'Woody Wetlands', 'Emergent Herbaceous Wetlands'};
manning = [0.038; 0.038; 0.040; 0.090; 0.120; 0.160; 0.027; 0.150; 0.120; 0.140; ...
           0.038; 0.115; 0.038; 0.038; 0.038; 0.035; 0.098; 0.068];


if exist('nlck1km.mat','file')
    load('nlck1km.mat');
else
    [nlcd,metadata] = readgeoraster('nlcd_2021_land_cover_l48_20230630.img');
    nlcd = nlcd(1:104412,1:161172);
    nlcd1km = NaN(104412/33,161172/33);
    for i = 1 : 3164
        disp(i);
        for j = 1 : 4884
            tmp = double(nlcd((i-1)*33+1+1:i*33,(j-1)*33+1:j*33));
            tmp(tmp == 0) = NaN;
            for k = 1 : length(vals)
                tmp(tmp == vals(k)) = manning(k);
            end
            nlcd1km(i,j) = nanmean(tmp(:));
        end
    end
end

w_bound_m = metadata.XWorldLimits(1);
e_bound_m = metadata.XWorldLimits(2);
s_bound_m = metadata.YWorldLimits(1);
n_bound_m = metadata.YWorldLimits(2);
x_step = metadata.CellExtentInWorldX;
y_step = metadata.CellExtentInWorldY;
% generate full list of x and y values
x_vals=w_bound_m+x_step/2:x_step:e_bound_m-x_step/2;
y_vals=n_bound_m-y_step/2:-y_step:s_bound_m+y_step/2;

x_vals = x_vals(1:161172);
y_vals = y_vals(1:104412);

x1km = NaN(size(nlcd1km,2),1);
y1km = NaN(size(nlcd1km,1),1);
for i = 1 : 3164
    y1km(i) = nanmean(y_vals((i-1)*33+1:i*33));
end
for j = 1 : 4884
    x1km(j) = nanmean(x_vals((j-1)*33+1:j*33));
end

[x1km,y1km] = meshgrid(x1km,y1km);

% project degree to m
% [x, y] = projfwd(metadata.ProjectedCRS,lat,lon);
