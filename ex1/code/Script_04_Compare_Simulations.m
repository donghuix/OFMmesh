clear;close all;clc;

[Output, exportfig] = SetupEnvironment();

smeshes =  {'Turning_30m', 'Turning_60m', 'Turning_90m', 'Turning_120m', ...
            'Turning_180m','Turning_240m','Turning_480m','Turning_990m'};

if exist('../data/ssims.mat','file')
    load('../data/ssims.mat');
else
    ssims = read_sims(Output,smeshes);
    save('../data/ssims.mat','ssims');
end

%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%
%-v-v-v-v-v-v-v-v-v-v-v-v-v- PLOT START HERE -v-v-v-v-v-v-v-v-v-v-v-v-v-%
%- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -%

tri    = ncread('../structure_meshes/Turning_30m.exo','connect1');
coordx = ncread('../structure_meshes/Turning_30m.exo','coordx');
coordy = ncread('../structure_meshes/Turning_30m.exo','coordy');
trix   = coordx(tri);
triy   = coordy(tri); 
load('MainChannel_poly.mat');
in = inpoly2([mean(trix)' mean(triy)'],[polyout.Vertices(:,1) polyout.Vertices(:,2)]);

letter = {'(a). ', '(b). ', '(c). ', '(d). ', '(e). ', '(f). ', '(g). ', '(h). ', '(i). '};
labels = {'30m','60m','90m','120m','180m','250m','500m','1,000m'};
load('colorblind_colormap.mat');

fig5 = figure(5); set(gcf,'Position',[10 10 1200 800]);
ax5(1) = subplot(3,4,[1 2 3 4]);
for i = 1 : length(smeshes)
    plot(ssims(i).t,ssims(i).q,'-','Color',colorblind(i,:),'LineWidth',2); hold on; grid on;
    set(gca,'FontSize',14); ylim([-100 8000]);
end
leg = legend(labels,'FontSize',13,'EdgeColor','none','color','none');
datetick('x','keeplimits');
ylabel('Discharge [m^{3}/s]','FontSize',15,'FontWeight','bold');

cmap = getPanoply_cMap('NEO_modis_lst');
%cmap = create_nonlinear_cmap(cmap,0,0.1,0.1,3);
for i = 1 : 8
    ax5(1+i) = subplot(3,4,4+i);
    patch(trix,triy,ssims(i).hmax_wholedomain,'LineStyle','none'); clim([0 12]);
    set(gca,'XTick',[],'YTick',[]); colormap(gca,cmap);
    if i == 8
        cb = colorbar('south');
    end
end
add_title(ax5(1),'(a). ',15,'in');

w  = ax5(1).Position(3);
mg = 0.01;
w0 = (w - mg*3)/4;

t = tight_axs(ax5,cb,ax5(1).Position(1),w0,mg,ssims,'hmax_wholedomain',letter,labels);
xlabel(cb,'Maximum inundation depth [m]','FontSize',15,'FontWeight','bold');

figure(501); set(gcf,'Position',[10 10 1200 800]);
for i = 1 : 8
    ax501(1+i) = subplot(3,4,4+i);
    patch(trix,triy,ssims(i).FI_wholedomain,'LineStyle','none'); clim([0 40]); hold on
    set(gca,'XTick',[],'YTick',[]); colormap(gca,cmap);
    if i == 8
        cb2 = colorbar('south');
    end
end
t2 = tight_axs(ax501,cb2,ax5(1).Position(1),w0,mg,ssims,'FI_wholedomain',letter,labels);
xlabel(cb2,'Maximum flood intensity [m^{2}/s]','FontSize',15,'FontWeight','bold');

fig502 = figure(502); set(gcf,'Position',[10 10 1200 800]);
for i = 1 : 8
    ax502(1+i) = subplot(3,4,4+i);
    patch(trix(:,in),triy(:,in),ssims(i).hmax_mainchannel,'LineStyle','none'); clim([0 12]); hold on
    set(gca,'XTick',[],'YTick',[]); colormap(gca,cmap);
    if i == 8
        cb3 = colorbar('south');
    end
end
t3 = tight_axs(ax502,cb3,ax5(1).Position(1),w0,mg,ssims,'hmax_mainchannel',letter,labels);
xlabel(cb3,'Maximum inundation depth [m]','FontSize',15,'FontWeight','bold');

if exportfig
    %exportgraphics(fig5,'Figure_5.pdf','ContentType','vector');
    exportgraphics(fig5,'Figure_5.jpg','Resolution',400);
end