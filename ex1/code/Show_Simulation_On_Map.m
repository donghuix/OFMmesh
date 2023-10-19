clear;close all;clc;

load('../data/depth.mat');
load('../data/ssims.mat');

proj = projcrs(32610);
[lat,lon] = projinv(proj,xbnd,ybnd);

[trilat,trilon] = projinv(proj,mesh.c_x,mesh.c_y);

cmap = getPanoply_cMap('NEO_omi_ozone_to3');
cmap(1:40,:) = [];
figure; set(gcf,'Position',[10 10 800 600]);
gx = geoaxes;
geoplot(gx,lat,lon,'r-','LineWidth',2);
geobasemap(gx,'satellite');

ax2 = axes;
ind = find(ssims(1).FI_wholedomain > 0.1);
patch(trilon(:,ind),trilat(:,ind),ssims(1).hmax_wholedomain(ind),'LineStyle','none','FaceAlpha',0.7); colormap(cmap); clim([0 12]);
xlim(gx.LongitudeLimits);
ylim(gx.LatitudeLimits);

% Set ax2 visibility to 'off'
ax2.Visible = 'off'; 
ax2.XTick = []; 
ax2.YTick = []; 

%exportgraphics(gcf,'Flood_on_map.jpg','Resolution',800);