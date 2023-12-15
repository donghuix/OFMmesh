clear;close all;clc;

addpath('/Users/xudo627/Developments/inpoly/');
addpath('/Users/xudo627/Developments/getPanoply_cMap/');
addpath('/Users/xudo627/Developments/petsc/share/petsc/matlab/');
addpath('/Users/xudo627/Developments/ofm_petsc/Matlab_Scripts/');

% Basins = {'Mackenzie', 'Mississippi', 'Orinoco', 'Amazon', 'Danube',  ...
%           'Volga', 'Ob', 'Godavari',  'Yangtze', 'Yenisey', ...
%           'Lena', 'Kolyma', 'Murray-Darling','Irrawaddy','Zaire'}; %'Pechora'
% ibasins = [5; 46; 149; 176; 43; 16; 253; 121; 97; 3; 252; 2; 228; 109; 156];
% basins  = shaperead('../major_basins_of_the_world_0_0_0/Major_Basins_of_the_World.shp');
% xb = basins(46).X';
% yb = basins(46).Y';
hybas_na_lev03 = shaperead('../hybas_na_lev01-12_v1c/hybas_na_lev03_v1c.shp');
xb = hybas_na_lev03(19).X';
yb = hybas_na_lev03(19).Y';

fname = '../data/hyd_na_dem_30s/hyd_na_dem_30s.tif';
dem = double(imread(fname));
dem(dem == 32767) = NaN;
I = geotiffinfo(fname); 
[x,y] =pixcenters(I);

% Simplify boundary
i = 1;
xb_simp = xb(1);
yb_simp = yb(1);
while i < length(xb)-1
    i = i + 1;
    dist =  sqrt((xb(i) - xb_simp(end))^2 + (yb(i) - yb_simp(end))^2)*100;
    if dist >= 1
        n = floor(dist);
        if n == 1
            xb_simp = [xb_simp; xb(i)];
            yb_simp = [yb_simp; yb(i)];
        elseif n > 1
            xb_interp = interp1([0; n],[xb_simp(end); xb(i)],[1:n]');
            yb_interp = interp1([0; n],[yb_simp(end); yb(i)],[1:n]');
            xb_simp = [xb_simp; xb_interp];
            yb_simp = [yb_simp; yb_interp];
        end
    end
end

i1 = min(find(y <= max(yb_simp)));
i2 = max(find(y >= min(yb_simp)));

j1 = min(find(x >= min(xb_simp)));
j2 = max(find(x <= max(xb_simp)));

x = x(j1-1:j2+1); 
y = y(i1-1:i2+1);
dem = dem(i1-1:i2+1,j1-1:j2+1);

% [m,n] = size(x);
% in = inpoly2([x(:) y(:)], [xb_simp yb_simp]);
% bcode = NaN(m,n);
% bcode(in) = 0;
% for i = 2 : m-1
%     disp(i);
%     for j = 2 : n-1
%         tmp = [bcode(i-1,j); bcode(i,j-1); bcode(i+1,j); bcode(i,j+1)];
%         if any(isnan(tmp))
%             bcode(i,j) = 1;
%         end
%     end
% end
% 
% x = x(in);
% y = y(in);
% bcode = bcode(in);
% 
% figure;
% plot(x(bcode == 1), y(bcode == 1),'bx','LineWidth',1); hold on;
% plot(xb_simp,yb_simp,'r-','LineWidth',2);

if exist('zinm.mat','file')
    load('zinm.mat');
else
    xm = mean([x(1:end-1)' x(2:end)'],2);
    ym = mean([y(1:end-1)' y(2:end)'],2);
    [xm,ym] = meshgrid(xm,ym);
    inm = inpoly2([xm(:) ym(:)], [xb yb]);
    zinm = griddata(x,y,dem,xm(inm),ym(inm));
    save('zinm.mat','xm','ym','inm','zinm');
end

[x,y] = meshgrid(x,y);
in = inpoly2([x(:) y(:)], [xb_simp yb_simp]);

bc_len = length(xb_simp);

if exist('zb_simp.mat','file')
    load('zb_simp.mat');
else
    zb_simp = griddata(x,y,dem,xb_simp,yb_simp);
    save('zb_simp.mat','zb_simp');
end
xc = [x(in);   xm(inm)]; 
yc = [y(in);   ym(inm)]; 
zc = [dem(in); zinm   ]; clear x y dem in;


polyout = polybuffer([xb_simp yb_simp],'lines',0.005);
ind = find(isnan(polyout.Vertices(:,1)));
figure;
plot(polyout.Vertices(1:ind-1,1),polyout.Vertices(1:ind-1,2),'k-'); hold on;
plot(polyout.Vertices(ind+1:end,1),polyout.Vertices(ind+1:end,2),'r-');

in_boundary = inpoly2([xc yc],polyout.Vertices(ind+1:end,:));

coordx = round([xb_simp; xc(in_boundary)],7);
coordy = round([yb_simp; yc(in_boundary)],7);
coordz = [zb_simp; zc(in_boundary)];

tmp = unique([coordx coordy],'rows','stable');
coordx = tmp(:,1); coordy = tmp(:,2); clear tmp;
bcode = zeros(length(coordx),1);
bcode(1:bc_len) = 1;

C = [ [1:bc_len]' [[2:bc_len]'; 1] ];
DT = delaunayTriangulation(coordx,coordy,C);
DT = DT.ConnectivityList;
in = inpoly2([nanmean(coordx(DT'))' nanmean(coordy(DT.'))'], [xb_simp yb_simp]);
DT(~in,:) = [];

area = polyarea(coordx(DT'),coordy(DT')).*1e4.*1e6; % about m^2
ismall = find(area <= 1e5);

figure;
plot(coordx(bc_len+1:end),coordy(bc_len+1:end),'k.'); hold on;
plot(coordx(1:bc_len),coordy(1:bc_len),'bo');
plot(nanmean(coordx(DT(ismall,:)')),nanmean(coordy(DT(ismall,:)')),'rx','LineWidth',2);
%plot([lon(DT(imin,:)'); lon(DT(imin,1)')],[lat(DT(imin,:)'); lat(DT(imin,1)')],'LineWidth',1.5);

cmap = getPanoply_cMap('GIST_earth');
inoutlet = inpoly2([nanmean(coordx(DT'))' nanmean(coordy(DT'))'], [[-89.4; -89.15; -89.15; -89.4; -89.4] [29.06; 29.06; 29.24; 29.24; 29.06]]);
figure;
patch(coordx(DT(inoutlet,:)'),coordy(DT(inoutlet,:)'),nanmean(zc(DT(inoutlet,:)'),1)); colorbar; hold on;
plot(nanmean(coordx(DT(inoutlet,:)')),nanmean(coordy(DT(inoutlet,:)')),'rx','LineWidth',1.5);
plot(xb,yb,'r--','LineWidth',1.5); 
xlim([-89.4 -89.15]); ylim([29.06 29.24]);
clim([0 10]);
% demcmap([0 200]);
colormap(cmap);

xout = -89.2886;
yout = +29.0925;
dist = (coordx(bcode == 1) - xout).^2 + (coordy(bcode == 1) - yout).^2;
iout = find(dist == min(dist));
figure;
plot(coordx(bcode == 1),coordy(bcode ==1),'bx'); hold on;
plot(coordx(iout), coordy(iout), 'ro', 'LineWidth',2);
bcode(iout) = 2;

b_tri = bcode(DT)';
PetscBinaryWrite('mississippi.bcode',b_tri(:));
write_exodus_file(coordx, coordy, coordz, DT, 'mississippi.exo');
% 
% if exist('mississippi_rivers.mat','file')
%     load('mississippi_rivers.mat');
% else
%     rivers = shaperead('../data/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na_shp/HydroRIVERS_v10_na.shp');
%     for i = 1 : length(rivers)
%         xr(i,1) = nanmean(rivers(i).BoundingBox(:,1));
%         yr(i,1) = nanmean(rivers(i).BoundingBox(:,2));
%     end
%     
%     river_in = inpoly2([xr yr],[xb yb]);
%     figure;
%     plot(xr(river_in),yr(river_in),'b.'); hold on;
%     rivers = rivers(river_in);
%     save('mississippi_rivers.mat','rivers');
% end
% 
% ORD_FLOW = [rivers(:).ORD_FLOW];
% cmap = getPanoply_cMap('GSFC_landsat_udf_density');
% cmap = flipud(cmap);
% 
% figure;
% plot(xb,yb,'r:','LineWidth',1); hold on;
% 
% ind = find(ORD_FLOW == 2);
% for i = 1 : length(ind)
%     plot(rivers(ind(i)).X,rivers(ind(i)).Y,'-','Color',[138,43,226]./255,'LineWidth',4); 
% end
% 
% ind = find(ORD_FLOW == 3);
% for i = 1 : length(ind)
%     plot(rivers(ind(i)).X,rivers(ind(i)).Y,'-','Color',[65,105,225]./255,'LineWidth',3); 
% end
% 
% ind = find(ORD_FLOW == 4);
% for i = 1 : length(ind)
%     plot(rivers(ind(i)).X,rivers(ind(i)).Y,'-','Color',[95,158,160]./255,'LineWidth',2); 
% end
% 
% ind = find(ORD_FLOW == 5);
% for i = 1 : length(ind)
%     plot(rivers(ind(i)).X,rivers(ind(i)).Y,'-','Color',[0,191,255]./255,'LineWidth',1); 
% end
% 
