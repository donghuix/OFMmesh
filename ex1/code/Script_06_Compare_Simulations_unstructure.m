clear;close all;clc;

[Output, exportfig] = SetupEnvironment();

umeshes =  {'Turning_30_60_1500_1000_fpx2', ...
            'Turning_30_60_1500_1000',      ...
            'Turning_30_90_1500_1000',      ...
            'Turning_60_90_1500_1000',      ...
            'Turning_90_90_1500_1000',      ...
            'Turning_30_90_1500_10000',     ...
            'Turning_30_90_1500_100000',    ... 
            'Turning_30_90_1500_1000000',   ...             
            'Turning_30_90_1500_NoRiver',   ...
            'Turning_River_Dam_Burn_1150'};

numc   = [664724, 427478, 254953, 197956, 185849, 148514, 109358, 58749, 14536, 2755];
if exist('../data/usims.mat','file')
    load('../data/usims.mat');
    load('../data/ssims.mat');
else
    usims = read_sims(Output,umeshes);
    save('../data/usims.mat','usims');
end

usims([4 5]) = [];
numc([4 5])  = [];

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%
%-v-v-v-v-v-v-v-v-v-v-v-v-v- PLOT START HERE -v-v-v-v-v-v-v-v-v-v-v-v-v-%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%

% batlow, scientific rainbow color scheme
load('/Users/xudo627/Developments/ofm_petsc/topotoolbox/colormaps/private/batlow/batlow.mat');

tri    = ncread('../structure_meshes/Turning_30m.exo','connect1');
coordx = ncread('../structure_meshes/Turning_30m.exo','coordx');
coordy = ncread('../structure_meshes/Turning_30m.exo','coordy');
trix   = coordx(tri);
triy   = coordy(tri); 
load('MainChannel_poly.mat');
in = inpoly2([mean(trix)' mean(triy)'],[polyout.Vertices(:,1) polyout.Vertices(:,2)]);

letter = {'(a). ', '(b). ', '(c). ', '(d). ', '(e). ', '(f). ', '(g). ', '(h). ', '(i). '};
labels = {'30m','Mesh 1','Mesh 2','Mesh 3','Mesh 4','Mesh 5','Mesh 6','Mesh 7','Mesh 8'};
load('colorblind_colormap.mat');
colorblind(3,:) = colorblind(12,:);
colorblind(9,:) = colorblind(11,:);

fig7 = figure(7); set(gcf,'Position',[10 10 1200 800]);
ax7(1) = subplot(3,4,[1 2 3 4]);
plot(ssims(1).t,ssims(1).q,'-','Color',colorblind(1,:),'LineWidth',3); hold on; grid on;
for i = 1 : length(usims)
    plot(usims(i).t,usims(i).q,'-','Color',colorblind(i+1,:),'LineWidth',2); 
    [qR2(i),qRMSE(i)] = estimate_evaluation_metric(ssims(1).q(:),usims(i).q(:));
    set(gca,'FontSize',14); ylim([-100 8000]);
end
leg = legend(labels,'FontSize',12,'EdgeColor','none','color','none');
datetick('x','keeplimits');
ylabel('Discharge [m^{3}/s]','FontSize',15,'FontWeight','bold');
cmap = getPanoply_cMap('NEO_modis_lst');

%cmap = create_nonlinear_cmap(cmap,0,0.1,0.1,3);
for i = 1 : 8
    ax7(1+i) = subplot(3,4,4+i);
    patch(trix,triy,usims(i).hmax_wholedomain,'LineStyle','none'); clim([0 12]);
    set(gca,'XTick',[],'YTick',[]); colormap(gca,cmap);
    if i == 8
        cb = colorbar('south');
    end
end
add_title(ax7(1),'(a). ',15,'in');

w  = ax7(1).Position(3);
mg = 0.01;
w0 = (w - mg*3)/4;

[t,R2,RMSE] = tight_axs(ax7,cb,ax7(1).Position(1),w0,mg,usims,'hmax_wholedomain',letter(2:end),labels(2:end),ssims);
xlabel(cb,'Maximum inundation depth [m]','FontSize',15,'FontWeight','bold');

fig702 = figure(702); set(gcf,'Position',[10 10 1200 800]);
for i = 1 : 8
    ax702(1+i) = subplot(3,4,4+i);
    patch(trix(:,in),triy(:,in),usims(i).hmax_mainchannel,'LineStyle','none'); clim([0 12]); hold on
    set(gca,'XTick',[],'YTick',[]); colormap(gca,cmap);
    if i == 8
        cb3 = colorbar('south');
    end
end
[t3,R23,RMSE3] = tight_axs(ax702,cb3,ax7(1).Position(1),w0,mg,usims,'hmax_mainchannel',letter,labels(2:end),ssims);
xlabel(cb3,'Maximum inundation depth [m]','FontSize',15,'FontWeight','bold');

fig8 = figure(8); set(gcf,'Position',[10 10 1200 400]);
numc2  = [750000,333333,187500,83333,49682,11718,2755];
structure = load('../data/structure_metrics.mat');
ax8(1) = subplot(1,3,1);
for i = 1 : 8
    semilogx(numc(i),qR2(i),'ks','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor',colorblind(i+1,:)); hold on; grid on;
end
semilogx(numc2,structure.qR2(2:end),'ko--','LineWidth',1.5); 
ylim([0 1]); set(gca,'FontSize',13);
add_title(ax8(1),'(a). Discharge');
ylabel('R^{2}','FontSize',15,'FontWeight','bold');

ax8(2) = subplot(1,3,2); 
for i = 1 : 8
    semilogx(numc(i),R2(i),'ks','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor',colorblind(i+1,:)); hold on; grid on;
end
semilogx(numc2,structure.R2(2:end),'ko--','LineWidth',1.5); 
ylim([0 1]); set(gca,'FontSize',13);
add_title(ax8(2),'(b). Maximum Inundation Depth');
xlabel('Number of cells','FontSize',15,'FontWeight','bold');

ax8(3) = subplot(1,3,3);
for i = 1 : 8
    h(i) = semilogx(numc(i),R23(i),'ks','LineWidth',1.5,'MarkerSize',12,'MarkerFaceColor',colorblind(i+1,:)); hold on; grid on;
end
h(length(numc)+1) = semilogx(numc2,structure.R23(2:end),'ko--','LineWidth',1.5); 
ylim([0 1]); set(gca,'FontSize',13);
add_title(ax8(3),'(c). Maximum Inundation Depth in Main Channel');

leg = legend(h,[labels(2:end), 'Uniform Meshes']);
leg.FontSize = 15;
ax8(2).Position(1) = ax8(2).Position(1) - 0.035;
ax8(3).Position(1) = ax8(3).Position(1) - 0.070;
leg.Position(1) = ax8(3).Position(1) + ax8(3).Position(3) + 0.01;
leg.Position(2) = ax8(3).Position(2);
leg.Position(4) = ax8(3).Position(4);

if exportfig
    %exportgraphics(fig5,'Figure_5.pdf','ContentType','vector');
    exportgraphics(fig7,'Figure_7.jpg','Resolution',400);
    exportgraphics(fig8,'Figure_8.pdf','ContentType','vector');
    exportgraphics(fig702,'Figure_702.jpg','Resolution',400);
end
