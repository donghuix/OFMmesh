clear;close all;clc;

[Output, exportfig] = SetupEnvironment();

labels = {'Less than Ankle','Ankle - Knee','Knee -Waist','Waist - Head','Above Head'};
% threshold for distance from a grid cell to outlet
thres  = 1000 : 1000 : 50000;

prefix = 'Turning_30m_SR';

if exist('../data/inundCF.mat','file')
    load('../data/inundCF.mat')
else
    [inundf,inundv,hdiff,inunds] = read_CF(Output, thres);
    save('../data/inundCF.mat','inunds','inundf','inundv','hdiff','-v7.3');
end
load('../data/depth.mat','xbnd','ybnd');

[x,y,z,b] = readbcode('../meshes/Turning_30m.bcode');
tri       = ncread('../meshes/Turning_30m.exo','connect1');

cmap = getPanoply_cMap('NEO_modis_lst');
cmap(1,:) = [1 1 1];
myblack  = [40 42 54]./255;
mygreen  = [80 250 123]./255;
mypurple = [149 128 255]./255;
%mypink   = [255 128 191]./255;
mypink   = [218 81 24]./255;

xmin = 3.210*1e6; 
xmax = 3.236*1e6;
ymin = 3.632*1e6;
ymax = 3.644*1e6;

load('../data/waterlevel.Outlet');
for i = 1 : size(waterlevel,1)
    twl(i,1) = datenum([waterlevel(i,1:5) 0]);
end
t0 = datenum(2017,8,26,0,0,0);
for i = 1 : 119
    tsim(i,1) = t0 + i*1/24;
end

fig4 = figure(4); set(gcf,'Position',[10 10 1200 900],'renderer','Painters')
ax(1) = subplot(2,2,1);
patch(x(tri),y(tri),hdiff,'LineStyle','none'); hold on; 
plot(xbnd,ybnd,'-','Color',mypurple,'LineWidth',2); hold on;
plot([xmin xmax xmax xmin xmin],[ymin ymin ymax ymax ymin],'--','Color',mypink,'LineWidth',3);
plot([3.22*1e6 3.23*1e6],[3.62*1e6 3.62*1e6], '-','Color',myblack,'LineWidth',3);
plot([3.22*1e6 3.22*1e6],[3.62*1e6 3.621*1e6],'-','Color',myblack,'LineWidth',3);
plot([3.23*1e6 3.23*1e6],[3.62*1e6 3.621*1e6],'-','Color',myblack,'LineWidth',3);
text(3.2225*1e6,3.622*1e6,'10km','Color',myblack,'FontSize',15,'FontWeight','bold');
% ylim([3.59*1e6 3.655*1e6])
% ax = axes('Position',[0.5 0.15 0.4 0.35]);
ax(2) = subplot(2,2,2);
patch(x(tri),y(tri),hdiff,'LineStyle','none'); hold on; 
plot(xbnd,ybnd,'-','Color',mypurple,'LineWidth',2); hold on;
clim([0 3]); colormap(cmap); cb = colorbar;
plot([3.229*1e6 3.234*1e6],[3.639*1e6 3.639*1e6], '-','Color',myblack,'LineWidth',3);
plot([3.229*1e6 3.229*1e6],[3.639*1e6 3.6393*1e6],'-','Color',myblack,'LineWidth',3);
plot([3.234*1e6 3.234*1e6],[3.639*1e6 3.6393*1e6],'-','Color',myblack,'LineWidth',3);
text(3.2305*1e6,3.6396*1e6,'5km','Color',myblack,'FontSize',15,'FontWeight','bold');
% Zoom in to river outlet
xlim([xmin xmax]); ylim([ymin ymax]);
set(gca,'XTick',[],'YTick',[]);
ylabel(cb,'\Delta h [m]','FontSize',16,'FontWeight','bold');

ax(3) = subplot(2,2,3);
dv = inundv(:,:,1) - inundv(:,:,2);
yyaxis right;
plot(twl,waterlevel(:,6),'-','Color',mypink,'LineWidth',2); grid on;
xlim([tsim(1) tsim(end)]);
datetick('x','keepticks');
set(gca,'FontSize',13);
ylabel('Water level [m]','FontSize',15,'FontWeight','bold');
yyaxis left;
plot(tsim,sum(dv,1),'b-','LineWidth',2);
ylabel('\Delta V [m^{3}/s]','FontSize',15,'FontWeight','bold');

for i = 1 : length(thres)
    dv(i,:) = sum(inundv(1:i,:,1) - inundv(1:i,:,2),1)./sum(inundv(1:i,:,2),1);
end

ax(4) = subplot(2,2,4);
plot(thres./1000,max(dv,[],2).*100,'k-','LineWidth',2); grid on;
set(gca,'FontSize',13);
ylabel('\Delta V [%]','FontSize',15,'FontWeight','bold');
xlabel('Distance to outlet [km]','FontSize',15,'FontWeight','bold');

ax(1).Position(1) = 0.050;
ax(1).Position(3) = 0.550;
ax(2).Position(1) = ax(1).Position(1)+ax(1).Position(3)+0.020;
ax(2).Position(3) = 0.275;
ax(3).Position(1) = ax(1).Position(1);
ax(3).Position(2) = ax(3).Position(2) + 0.05;
w = (ax(2).Position(1) + ax(2).Position(3) - ax(1).Position(1) - 0.1) / 2;
ax(3).Position(3) = w;
ax(4).Position(1) = ax(3).Position(1) + w + 0.1;
ax(4).Position(3) = w;
ax(4).Position(2) = ax(4).Position(2) + 0.05;

Labels = {'(a)','(b)','(c)','(d)'};
for i = 1 : 4
    dim = [ax(i).Position(1) ax(i).Position(2)+ax(i).Position(4)-0.1 .1 .1];
    annotation('textbox',dim,'String',Labels{i},'FitBoxToText','on',    ...
               'FontSize',18,'FontWeight','bold','EdgeColor','none');
end

if exportfig
    exportgraphics(fig4,'Figure_4.jpg','Resolution',400);
end