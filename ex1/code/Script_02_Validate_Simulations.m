clear;close all;clc;

[Output, proj] = SetupEnvironment();
names = {'Turning_30m','Turning_30m_SR','Turning_30m_noBC'};
if exist('../data/depth.mat','file')
    load('../data/depth.mat');
else
    [xbnd,ybnd,gages,hwm,mesh,dams] = read_data(names,Output,proj,200);
    save('../data/depth.mat','xbnd','ybnd','gages','hwm','mesh','dams');
end

Pix_SS = get(0,'screensize');
if Pix_SS(3)/Pix_SS(4) > 1.6
    fw(1) = 1200;
    fd(1) = 800;
    fw(2) = 1800;
    fd(2) = 1000;
end

figure(1); set(gcf,'Position',[10 10 1200 800]);
cmap = getPanoply_cMap('GIST_earth');
patch(mesh.c_x,mesh.c_y,mesh.z_tri,'LineStyle','none'); hold on;
colormap(cmap); 
cb = colorbar;
cb.FontSize = 15; 

plot([hwm.X],[hwm.Y],'wx','LineWidth',1);
set(gca,'Color',[0.75 0.75 0.75]);

plot(dams(1).X,dams(1).Y,'k-','LineWidth',2); 
plot(dams(2).X,dams(2).Y,'k-','LineWidth',2); 


figure(2); set(gcf,'Position',[10 10 1200 800]);
k = 1;
plot(xbnd,ybnd,'k-','LineWidth',2); hold on; grid on;
for i = 1:length(gages)
    if ~isempty(gages(i).wl)
        figure(2);
        axs(k) = subplot(4,5,k);
        plot(gages(i).tsim,gages(i).wl,'k-','LineWidth',2); hold on; grid on;

        %plot(tsim,hsim(ind1,:),'r-','LineWidth',2); hold on;
        plot(gages(i).tsim,gages(i).sim(1).h,'b-','LineWidth',2); 
        plot(gages(i).tsim,gages(i).sim(2).h,'r--','LineWidth',2); 
        plot(gages(i).tsim,gages(i).sim(3).h,'g--','LineWidth',2); 
        if k == 19
            leg = legend('USGS OBS','Sim','Sim with uniform rainfall','Sim with critical BC');
            leg.FontSize   = 15;
            leg.FontWeight = 'bold';
        end
        set(gca,'FontSize',12);
        datetick('x','dd/mm','keeplimits');
        add_title(gca,['Gauge#' num2str(k)]);
        figure(1);
        plot(gages(i).X,gages(i).Y,'*','Color',[255,69,0]./255,'LineWidth',1.5);
        t(k) = text(gages(i).X+500,gages(i).Y,['#' num2str(k)],'Color', ...
                    [255,69,0]./255,'FontSize',14,'FontWeight','bold');
        k = k + 1;
    end
end
leg.Position(1) = axs(15).Position(1);
leg.Position(2) = axs(19).Position(2) + axs(19).Position(4) - leg.Position(4);


hmax = 3; hmin = 0.1;
figure(3);
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