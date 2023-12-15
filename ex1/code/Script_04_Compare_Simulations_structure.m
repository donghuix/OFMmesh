clear;close all;clc;

[Output, exportfig] = SetupEnvironment();

smeshes =  {'Turning_30m', 'Turning_60m', 'Turning_90m', 'Turning_120m', ...
            'Turning_180m','Turning_240m','Turning_480m','Turning_990m'};

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

if exist('../data/ssims.mat','file')
    load('../data/ssims.mat');
else
    ssims = read_sims(Output,smeshes);
    save('../data/ssims.mat','ssims');
end

if exist('../data/usims.mat','file')
    load('../data/usims.mat');
else
    usims = read_sims(Output,umeshes);
    save('../data/usims.mat','usims');
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

letter = {'(b). ', '(c). ', '(d). ', '(e). ', '(f). ', '(g). ', '(h). ', '(i). '};
labels = {'30m','60m','90m','120m','180m','250m','500m','1,000m'};
load('colorblind_colormap.mat');

fig5 = figure(5); set(gcf,'Position',[10 10 1200 800]);
ax5(1) = subplot(3,4,[1 2 3 4]);
for i = 1 : length(smeshes)
    plot(ssims(i).t,ssims(i).q,'-','Color',colorblind(i,:),'LineWidth',2); hold on; grid on;
    if i > 1
        [qR2(i),qRMSE(i)] = estimate_evaluation_metric(ssims(1).q(:),ssims(i).q(:));
    end
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

[t,R2,RMSE] = tight_axs(ax5,cb,ax5(1).Position(1),w0,mg,ssims,'hmax_wholedomain',letter,labels);
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
[t3,R23,RMSE3] = tight_axs(ax502,cb3,ax5(1).Position(1),w0,mg,ssims,'hmax_mainchannel',letter,labels);
xlabel(cb3,'Maximum inundation depth [m]','FontSize',15,'FontWeight','bold');

if exportfig
    %exportgraphics(fig5,'Figure_5.pdf','ContentType','vector');
    exportgraphics(fig5,'Figure_5.jpg','Resolution',400);
end

save('../data/structure_metrics.mat','qR2','qRMSE','R2','RMSE','R23','RMSE3');

% function t = tight_axs(ax,cb,x0,w0,mg,ssims,varname,letter,labels)
%     %num = {'~3,000,000','750,000','333,333','187,500','83,333','49,682','11,718','2,755'};
%     for i = 1 : 4
%         % second row
%         ax(1+i).Position(3) = w0;
%         ax(1+i).Position(1) = x0 + (i-1)*(w0+mg);
%         ax(1+i).Position(2) = ax(1+i).Position(2) + 0.04;
%         % third row
%         ax(5+i).Position(3) = w0;
%         ax(5+i).Position(1) = x0 + (i-1)*(w0+mg);
%         ax(5+i).Position(2) = ax(1+i).Position(2) - mg - ax(5+i).Position(4);
%     end
%     for i = 1 : 8
%         add_title(ax(i+1),[letter{i} labels{i}],15,'in');
%         if i > 1
%             [R2,RMSE,~,PBIAS] = estimate_evaluation_metric(ssims(1).(varname),ssims(i).(varname));
%             R2    = round(R2,2);
%             RMSE  = round(RMSE,2);
%             str   = {['R^{2} = ' num2str(R2)], ['RMSE= ' num2str(RMSE)]};
%             t(i)  = add_annot(ax(i+1),str,12,'lr');
%         end
% %         tt(i) = add_annot(ax(i+1),num{i},12,'ur');
% %         tt(i).Color = 'r';
%     end
% 
%     cb.Position(1) = ax(6).Position(1);
%     cb.Position(2) = ax(6).Position(2) - 0.075;
%     cb.Position(3) = ax(9).Position(1) + ax(9).Position(3) - ax(6).Position(1);
%     cb.Position(4) = 0.02;
%     cb.FontSize = 14;
% end