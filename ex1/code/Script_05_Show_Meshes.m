clear;close all;clc;

[Output, exportfig] = SetupEnvironment();

meshes =  {'Turning_30_60_1500_1000_fpx2', ...
          'Turning_30_60_1500_1000',      ...
          'Turning_30_90_1500_1000',      ...
          'Turning_60_90_1500_1000',      ...
          'Turning_90_90_1500_1000',      ...
          'Turning_30_90_1500_10000',     ...
          'Turning_30_90_1500_100000',    ... 
          'Turning_30_90_1500_1000000',   ...             
          'Turning_30_90_1500_NoRiver',   ...
          'Turning_River_Dam_Burn_1150'};

meshes([4 5]) = [];
letter = {'(a). ', '(b). ', '(c). ', '(d). ', '(e). ', '(f). ', '(g). ', '(h). ', '(i). '};
labels = {'Mesh 1','Mesh 2','Mesh 3','Mesh 4','Mesh 5','Mesh 6','Mesh 7','Mesh 8'};
num    = {'664,724','427,478','254,953','148,514','109,358','58,749','14,536','~3,000'};

% Projection coordinates
proj = projcrs(32610);
% Read processed DEM
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

DEM = GRIDobj(xdem,ydem,dem);
FD  = FLOWobj(DEM);
ACC = flowacc(FD);
xpos = xdem(1,:)';
ypos = flipud(ydem(:,1))';
[XPOS,YPOS] = meshgrid(xpos,ypos); XPOS = XPOS'; YPOS = YPOS';
[m,n] = size(xdem);
in    = inpoly2([xdem(:) ydem(:)],[xbnd ybnd]);

f = flowacc(FD);
f.Z(~in) = 1;

%cmap = [4 89 147; 219 96 0; 180 12 13; 219 219 23]./255;

cmap = viridis();
fig6 = figure(6); set(gcf,'Position',[10 10 1200 800]);
for i = 1 : length(meshes)
    ax6(i) = subplot(3,4,i);
    mesh = meshes{i};
    [x,y,z,b] = readbcode(['../meshes/' mesh '.bcode']);
    tri = ncread(['../meshes/' mesh '.exo'],'connect1');
    xtri = x(tri);
    ytri = y(tri);
    area = polyarea(xtri,ytri);
    
    patch(xtri,ytri,sqrt(area),'LineStyle','none');
    set(gca,'ColorScale','log');
    set(gca,'XTick',[],'YTick',[]);
    clim([30 700]); colormap(gca,cmap);
    if i == 8
        cb = colorbar('south');
    end
end

w  = ax6(4).Position(1) + ax6(4).Position(3) - ax6(1).Position(1);
mg = 0.01;
w0 = (w - mg*3)/4;
x0 = ax6(1).Position(1);

for i = 1 : 4
    % second row
    ax6(i).Position(3) = w0;
    ax6(i).Position(1) = x0 + (i-1)*(w0+mg);
    ax6(i).Position(2) = ax6(i).Position(2) + 0.04;
    % third row
    ax6(4+i).Position(3) = w0;
    ax6(4+i).Position(1) = x0 + (i-1)*(w0+mg);
    ax6(4+i).Position(2) = ax6(i).Position(2) - mg - ax6(4+i).Position(4);

    add_title(ax6(i),  [letter{i} labels{i}],15,'in');
    add_title(ax6(4+i),[letter{4+i} labels{4+i}],15,'in');
    tt(i) = add_annot(ax6(i),num{i},12,'ur');
    tt(4+i) = add_annot(ax6(4+i),num{4+i},12,'ur');
    tt(i).Color = 'r';
    tt(4+i).Color = 'r';
end

cb.Position(1) = ax6(5).Position(1);
cb.Position(2) = ax6(5).Position(2) - 0.075;
cb.Position(3) = ax6(8).Position(1) + ax6(8).Position(3) - ax6(5).Position(1);
cb.Position(4) = 0.02;
cb.FontSize = 14;
cb.XTick = [30 90 700];
cb.XTickLabel = {'30', '90', '700'};
xlabel(cb,'Cell resolution [m]','FontSize',15,'FontWeight','bold');

if exportfig
    %exportgraphics(fig5,'Figure_5.pdf','ContentType','vector');
    exportgraphics(fig6,'Figure_6.jpg','Resolution',400);
end
% for im = 5
%     strs = strsplit(meshes{im},'_');
%      % parameters to generate mesh
%     hmain = str2double(strs{2});    % Main channle cell density
%     htrib = str2double(strs{3});    % Tributary cell density 
%     hhill = str2double(strs{4});    % Hillslope cell density
%     raccu = str2double(strs{5});    % Accumulation threshold for river
%     fp_width1  = 250;               % main channel floodplain width  [m]
%     fp_width2  = 150;               % tributries floodplain width    [m]
%     rm_width   = 180;               % river mouth width [m]  
%     dhdx       = 1;                 % Marche gradient smooth
%     add_dam    = 0;
% 
%     if length(strs) == 6
%         fp_width1 = fp_width1*2;
%         fp_width2 = fp_width2*2;
%     end
%     
%     %subplot(5,3,i);
%     
%     hfun  = ones(m,n).*0;
%     
%     smobj = STREAMobj(FD,f>=raccu);
%     order = smobj.orderednanlist;
%     streams_seg = [];
%     idxs  = find(isnan(order));
%     acc   = NaN(length(idxs),1);
%     for i = 1 : length(idxs)
%         if i == 1
%             acc(i) = max(f.Z(smobj.IXgrid(order(1:idxs(i)-1))));
%         else
%             acc(i) = max(f.Z(smobj.IXgrid(order(idxs(i-1)+1:idxs(i)-1))));
%         end
%     end
%     imain = find(acc == max(acc));
%     for i = 1 : length(idxs)
%         if i == 1
%             smx = smobj.x(order(1:idxs(i)-1));
%             smy = smobj.y(order(1:idxs(i)-1));
%             
%         else
%             smx = smobj.x(order(idxs(i-1)+1:idxs(i)-1));
%             smy = smobj.y(order(idxs(i-1)+1:idxs(i)-1));
%             acc = f.Z(smobj.IXgrid(order(idxs(i-1)+1:idxs(i)-1)));
%         end
%         polyout = polybuffer([smx,smy],'lines',fp_width2);
%         in = inpoly2([xdem(:) ydem(:)],polyout.Vertices);
%         hfun(in) = 1;
%         if i == imain
%             irm = find(smx < 3.195*1e6);
%             smx(irm) = [];
%             smy(irm) = [];
%             polyout = polybuffer([smx,smy],'lines',fp_width1);
%             main_in = inpoly2([xdem(:) ydem(:)],polyout.Vertices);
%         end
%         %plot(polyout); hold on;
%     end
%     hfun(main_in) = 2;
% 
%     for i = 1 : length(smobj.ix)
%         x1 = smobj.x(smobj.ix(i));
%         x2 = smobj.x(smobj.ixc(i));
%         y1 = smobj.y(smobj.ix(i));
%         y2 = smobj.y(smobj.ixc(i));
%         %plot([x1 x2], [y1 y2],'b.-','LineWidth',1); hold on;
%     end
%     %plot(xbnd,ybnd,'g.-','LineWidth',2); grid on;
% 
%     d1 = polybuffer([dam1.X' dam1.Y'],'lines',50);
%     d2 = polybuffer([dam2.X' dam2.Y'],'lines',50);
%     in1 = inpoly2([xdem(:) ydem(:)],d1.Vertices);
%     in2 = inpoly2([xdem(:) ydem(:)],d2.Vertices);
%     hfun(in1) = 3;
%     hfun(in2) = 3;
%     
%     in = inpoly2([xdem(:) ydem(:)],[xbnd ybnd]);
%     hfun(~in) = NaN;
%     hfun = flipud(hfun);
%     
%     imAlpha = ones(size(hfun));
%     imAlpha(isnan(hfun)) = 0;
%     imagesc([XPOS(1,1), XPOS(end,end)],[YPOS(1,1) YPOS(end,end)],hfun,'AlphaData',imAlpha);
%     set(gca,'YDir','normal');
%     colormap(cmap);
%     clim([-0.5 3.5]); colormap;
% end

