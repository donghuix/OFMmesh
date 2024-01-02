clear;close all;clc;

[Output, exportfig] = SetupEnvironment();
proj = projcrs(32610);

if exist('aux.mat','file')
    load('aux.mat');
else
    [xbnd,ybnd,gages,hwm,dams] = read_data();
    save('aux.mat','xbnd','ybnd','gages','hwm','dams');
end
names = {'Turning_30m','Turning_30m_SR','Turning_30m_noBC_SR'};
if exist('../data/depth.mat','file')
    load('../data/depth.mat');
else
    [gages,hwm,mesh] = read_depth(gages,hwm,names,Output,100);
    save('../data/depth.mat','xbnd','ybnd','gages','hwm','mesh','dams');
end

Pix_SS = get(0,'screensize');
if Pix_SS(3)/Pix_SS(4) > 1.6
    fw(1) = 1200;
    fd(1) = 800;
    fw(2) = 1800;
    fd(2) = 1000;
end

fig1 = figure(1); set(gcf,'Position',[10 10 1200 600]);
cmap = getPanoply_cMap('GIST_earth');
[lat,lon] = projinv(proj,mesh.c_x,mesh.c_y);
%patch(mesh.c_x,mesh.c_y,mesh.z_tri,'LineStyle','none'); hold on;
patch(lon,lat,mesh.z_tri,'LineStyle','none'); hold on;
set(gca,'Color',[0.75 0.75 0.75],'FontSize',15);
colormap(cmap); 
cb = colorbar;
cb.FontSize = 15; 
xlabel(cb,'Elevation [m]','FontSize',18,'FontWeight','bold');
% Add locations of HWM
[lat,lon] = projinv(proj,[hwm.X],[hwm.Y]);
plot(lon,lat,'wx','LineWidth',1.5,'MarkerSize',8);
% Add dams
[lat,lon] = projinv(proj,dams(1).X,dams(1).Y);
plot(lon,lat,'k-','LineWidth',2); 
[lat,lon] = projinv(proj,dams(2).X,dams(2).Y);
plot(lon,lat,'k-','LineWidth',2); 
ylim([29.68 29.96]); xlim([-95.975 -95.25]);

fig2 = figure(2); set(gcf,'Position',[10 10 1200 800],'renderer','Painters');
k = 1;
plot(xbnd,ybnd,'k-','LineWidth',2); hold on; grid on;
for i = 1:length(gages)
    if ~isempty(gages(i).wl)
        figure(2);
        axs(k) = subplot(4,5,k);
        plot(gages(i).tsim,gages(i).wl,'k-','LineWidth',2); hold on; grid on;

        %plot(tsim,hsim(ind1,:),'r-','LineWidth',2); hold on;
        plot(gages(i).tsim,gages(i).sim(1).h,'b--','LineWidth',2); 
%         plot(gages(i).tsim,gages(i).sim(2).h,'r--','LineWidth',2); 
%         plot(gages(i).tsim,gages(i).sim(3).h,'g--','LineWidth',2); 
        [R2(k),RMSE(k),NSE(k)] = estimate_evaluation_metric(gages(i).wl(:),gages(i).sim(1).h(:));
        if k == 19
            %leg = legend('USGS OBS','Sim','Sim with uniform rainfall','Sim with critical BC');
            leg = legend('Observation','Simulation');
            leg.FontSize   = 15;
            leg.FontWeight = 'bold';
        end
        datetick('x','dd'); xlim([gages(i).tsim(1) gages(i).tsim(end)]);
        set(gca,'FontSize',12);
        add_title(gca,['Gauge#' num2str(k)]);
        strs = {['\rho = ' num2str(round(sqrt(R2(k)),2)) ','], ['NSE = ' num2str(NSE(k))]};
        t2(k) = add_title(axs(k),strs,13,'in');
        %t2(k).Position(1) = t2(k).Position(1) -0.5 + t2(k).Position(3)/2;
        pos = get(axs(k),'Position');
        t2(k).Position(1) = t2(k).Position(1) + 0.03;
        t2(k).Position(2) = pos(2)-0.04;
        t2(k).Color = 'b'; 

        figure(1);
        [lat,lon] = projinv(proj,gages(i).X,gages(i).Y);
        plot(lon,lat,'*','Color',[255,69,0]./255,'LineWidth',1.5);
        [lat,lon] = projinv(proj,gages(i).X+500,gages(i).Y);
        t(k) = text(lon,lat,['#' num2str(k)],'Color', ...
                    [255,69,0]./255,'FontSize',14,'FontWeight','bold');
        
        k = k + 1;
    end
end
leg.Position(1) = axs(15).Position(1);
leg.Position(2) = axs(19).Position(2) + axs(19).Position(4) - leg.Position(4);

han2=axes(fig2,'visible','off'); 
han2.YLabel.Visible='on';
ylabel(han2,'Water level [m]','FontSize',15,'FontWeight','bold');
han2.Position(1) = han2.Position(1) - 0.01;

hmax = 3; hmin = 0.1;
fig3 = figure(3);
loglog([hmin hmax],[hmin hmax],'-','Color',[220,20,60]./255,'LineWidth',3);hold on; grid on;
r = 1/2;
loglog([hmin/r hmax],[hmin hmax*r],'--','Color',[220,20,60]./255,'LineWidth',2);
loglog([hmin hmax*r],[hmin/r hmax],'--','Color',[220,20,60]./255,'LineWidth',2);
r = 3/4;
loglog([hmin/r hmax],[hmin hmax*r],':','Color',[220,20,60]./255,'LineWidth',2);
loglog([hmin hmax*r],[hmin/r hmax],':','Color',[220,20,60]./255,'LineWidth',2);
loglog([hwm.h],[hwm.sim(1).h],'ko','MarkerSize',10,'MarkerFaceColor',[00,149,237]./255,'LineWidth',2); 


xlim([hmin hmax]);
ylim([hmin hmax]);
set(gca,'XTick',[0.1:0.1:0.5 1 2 3],'YTick',[0.1:0.1:0.5 1 2 3],'FontSize',14);
set(gca,'FontSize',13);
xlabel('High Water Marks [m]','FontSize',15,'FontWeight','bold');
ylabel('Simulation [m]','FontSize',15,'FontWeight','bold');

if exportfig
    % exportgraphics(fig1,'Figure_1.pdf','ContentType','image');
    % To crop the margins, use the following command
    % pdfcrop --margins '-85 -35 -75 -35' Figure_1.pdf
    exportgraphics(fig1,'Figure_1.jpg','Resolution',400);
    exportgraphics(fig2,'Figure_3.pdf','ContentType','vector');
    exportgraphics(fig3,'Figure_4.pdf','ContentType','vector');
end

 % Projection coordinates
figure;
[lat,lon] = projinv(proj,xbnd,ybnd);
geoplot(lat,lon,'b-','LineWidth',2); hold on;
geobasemap topographic;

x05 = -96 + 0.5/2 : 0.5 : -95 - 0.5/2;
y05 = 29 + 0.5/2 : 0.5 : 30 - 0.5/2;
[x05,y05] = meshgrid(x05,y05);
[xv,yv] = xc2xv(x05(:),y05(:),0.5,0.5);

for i = 2
    geoplot([yv(:,i); yv(1,i)],[xv(:,i); xv(1,i)],'r-','LineWidth',1);
end

Rfiles = dir('/Users/xudo627/Developments/ofm_petsc/Input/Rain/Turning/*.asc');
R      = NaN(38,77,length(Rfiles));
for i = 1 : length(Rfiles)
    Rfilename = fullfile(Rfiles(i).folder,Rfiles(i).name);
    [R(:,:,i), X, Y] = ascread(Rfilename,[]);
end
% figure;
% for i = 1 : length(Rfiles)
%     subplot(4,6,i);
%     imagesc(R(:,:,i));
% end

Rtri = griddata(X,Y,mean(R,3),mesh.x_tri',mesh.y_tri');
figure; set(gcf,'Position',[10 10 800 400])
title('2017-Aug-27 mean rainfall rate [mm/hour]','FontSize',15,'FontWeight','bold');
patch(mesh.c_x,mesh.c_y,Rtri,'LineStyle','none'); cb = colorbar; axis equal;
cb.FontSize=13;
set(gca,'XTick',[],'YTick',[]);
exportgraphics(gcf,'Figure_y.jpg','Resolution',400);

data = dlmread('Turning.Rain','\t',1, 0);
t = datenum(data(:,1),data(:,2),data(:,3),data(:,4),ones(168,1),ones(168,1));
figure;
bar(t,data(:,5),'k');
datetick('x','mm/dd','keepticks'); grid on;
xlim([t(1) t(end-36)]);
set(gca,'FontSize',13);
ylabel('Rainfall rate [mm/hour]','FontSize',15,'FontWeight','bold');
exportgraphics(gcf,'Figure_z.pdf','ContentType','vector');
% mdl = fitlm([hwm.h],[hwm.sim(1).h]);
% a = mdl.Coefficients.Estimate(1);
% b = mdl.Coefficients.Estimate(2);
% plot([hmin hmax],[hmin*b+a hmax*b+a],'-','Color',[220,20,60]./255,'LineWidth',2);
% ci95 = coefCI(mdl);
% 
% y1 = hmin*ci95(2,1) + ci95(1,1);
% y2 = hmax*ci95(2,1) + ci95(1,1);
% y3 = hmax*ci95(2,2) + ci95(1,2);
% y4 = hmin*ci95(2,2) + ci95(1,2);
% fill([hmin hmax hmax hmin],[y1 y2 y3 y4],[220,20,60]./255,'FaceAlpha',0.2,'EdgeColor','none');

% ind1 = find([hwm.sim(1).h] >= hwm(1).h.*ci95(2,1) + ci95(1,1) & ...
%             [hwm.sim(1).h] <= hwm(1).h.*ci95(2,2) + ci95(1,2));

% figure;
% load('../smobj.mat');
% plot(xbnd,ybnd,'k-','LineWidth',2); hold on; grid on;
% for i = length(sites)
%     [x,y] = projfwd(proj,lat(i),lon(i));
%     scatter(x,y,72,'s','MarkerFaceColor','r','MarkerEdgeColor','k','LineWidth',2);
%     plot(x_tri(idxusgs(i,:)),y_tri(idxusgs(i,:)),'r.')
% end
% for i = 1 : length(smobj.ix)
%     x1 = smobj.x(smobj.ix(i));
%     x2 = smobj.x(smobj.ixc(i));
%     y1 = smobj.y(smobj.ix(i));
%     y2 = smobj.y(smobj.ixc(i));
% 
%     plot([x1 x2], [y1 y2],'b.-','LineWidth',1); hold on;
% end
% 
% figure;
% patch(c_x,c_y,z_tri,'LineStyle','none'); colorbar; hold on;
% plot(xbnd,ybnd,'g-','LineWidth',2); hold on; grid on;
% for i = length(sites)
%     [x,y] = projfwd(proj,lat(i),lon(i));
%     plot(x,y,'rs','LineWidth',2);
% end