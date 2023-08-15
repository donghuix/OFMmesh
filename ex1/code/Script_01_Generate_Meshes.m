clear;close all;clc;

% NOTES:
% (1). Burn in streamline geometries will result in tiny triangles, e.g., area < 0.1 m^2.
% (2). Burn in dam geometries doesn't improve the simulation, need to zoom in dam regions with high resolution.
% (3). 

% parameters to generate mesh
river_accs = [1e3];   % [2.5e3 1e4 2.5e4 5e4 1e5 2.5e5 3e5 4e5 5e5];
fp_width1  = 250;     % main channel floodplain width  [m]
fp_width2  = 150;     % tributries floodplain width    [m]
rm_width   = 180;     % river mouth width [m]
hmain      = 30;      % Main channle cell density
htrib      = 90;      % Tributary cell density 
hhill      = 1500;    % Hillslope cell density
dhdx       = 1;       % Marche gradient smooth
add_dam    = 0;

if hmain > 1000
    add_stream = 1;
else
    add_stream = 0;
end

[Output, proj] = SetupEnvironment();

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


figure;
imagesc([xdem(1,1) xdem(end,end)],[ydem(1,1) ydem(end,end)],dem);
set(gca,'YDir','normal'); colormap(gca,jet); hold on; 
plot(dam1.X,dam1.Y,'g-','LineWidth',1);
plot(dam2.X,dam2.Y,'g-','LineWidth',1); caxis([0 30]);

DEM = GRIDobj(xdem,ydem,dem);
FD  = FLOWobj(DEM);
ACC = flowacc(FD);

xpos = xdem(1,:)';
ypos = flipud(ydem(:,1))';
[XPOS,YPOS] = meshgrid(xpos,ypos); XPOS = XPOS'; YPOS = YPOS';
[m,n] = size(xdem);
hfun = ones(m,n).*hhill;
in = inpoly2([xdem(:) ydem(:)],[xbnd ybnd]);

k = 1;
figure;
for river_acc = river_accs
    subplot(1,2,k);
    f = flowacc(FD);
    f.Z(~in) = 1;
    smobj = STREAMobj(FD,f>=river_acc);
    order = smobj.orderednanlist;
    streams_seg = [];
    idxs  = find(isnan(order));
    acc   = NaN(length(idxs),1);
    for i = 1 : length(idxs)
        if i == 1
            acc(i) = max(f.Z(smobj.IXgrid(order(1:idxs(i)-1))));
        else
            acc(i) = max(f.Z(smobj.IXgrid(order(idxs(i-1)+1:idxs(i)-1))));
        end
    end
    imain = find(acc == max(acc));

    for i = 1 : length(idxs)
        if i == 1
            smx = smobj.x(order(1:idxs(i)-1));
            smy = smobj.y(order(1:idxs(i)-1));
            
        else
            smx = smobj.x(order(idxs(i-1)+1:idxs(i)-1));
            smy = smobj.y(order(idxs(i-1)+1:idxs(i)-1));
            acc = f.Z(smobj.IXgrid(order(idxs(i-1)+1:idxs(i)-1)));
        end
        polyout = polybuffer([smx,smy],'lines',fp_width2);
        in = inpoly2([xdem(:) ydem(:)],polyout.Vertices);
        hfun(in) = htrib;
        if i == imain
            irm = find(smx < 3.195*1e6);
            smx(irm) = [];
            smy(irm) = [];
            polyout = polybuffer([smx,smy],'lines',fp_width1);
            main_in = inpoly2([xdem(:) ydem(:)],polyout.Vertices);
        end
        plot(polyout); hold on;
    end
    hfun(main_in) = hmain;

    for i = 1 : length(smobj.ix)
        x1 = smobj.x(smobj.ix(i));
        x2 = smobj.x(smobj.ixc(i));
        y1 = smobj.y(smobj.ix(i));
        y2 = smobj.y(smobj.ixc(i));

        plot([x1 x2], [y1 y2],'b.-','LineWidth',1); hold on;
    end
    plot(xbnd,ybnd,'g.-','LineWidth',2); grid on;
    k = k + 1;
end

d1 = polybuffer([dam1.X' dam1.Y'],'lines',50);
d2 = polybuffer([dam2.X' dam2.Y'],'lines',50);
in1 = inpoly2([xdem(:) ydem(:)],d1.Vertices);
in2 = inpoly2([xdem(:) ydem(:)],d2.Vertices);
hfun(in1) = 30;
hfun(in2) = 30;

in = inpoly2([xdem(:) ydem(:)],[xbnd ybnd]);
hfun(~in) = hhill;
hfun = flipud(hfun);

rootpath = '../meshes/' ;
name = ['Turning_' num2str(hmain) '_' num2str(htrib) '_' num2str(hhill) '_'  num2str(river_acc)];

opts.geom_file = fullfile(rootpath,[name,'.msh']) ;
opts.jcfg_file = fullfile(rootpath,[name,'.jig']) ;    
opts.mesh_file = fullfile(rootpath,[name,'.msh']) ;
opts.hfun_file = fullfile(rootpath,[name,'-hfun.msh']);
initjig ;                           % init jigsaw

% Add geometry for boundary (bd), streamline (sm), two reservoirs (d1, d2)
bd  = [xbnd(1:end-1) ybnd(1:end-1) ones(length(xbnd)-1,1)];
sm  = [smobj.x smobj.y zeros(length(smobj.x),1)];
d1  = [dam1.X(1:end-1)' dam1.Y(1:end-1)' zeros(length(dam1.X)-1,1)];
d2  = [dam2.X(1:end-1)' dam2.Y(1:end-1)' zeros(length(dam2.X)-1,1)];

% Make sure the geometry is inside the watershed boundary
[xi,yi] = polyxpoly(d2(:,1),d2(:,2),bd(:,1),bd(:,2));
circ = polybuffer([xi yi],'points',100);
[xi,yi] = polyxpoly(d2(:,1),d2(:,2),circ.Vertices(:,1),circ.Vertices(:,2));
d2(end,1) = xi;
d2(end,2) = yi;

irm = [];
for i = 1 : length(smobj.ix)
    disp(['i = ' num2str(i) '/' num2str(length(smobj.ix))]);
    smx = smobj.x([smobj.ix(i) smobj.ixc(i)]);
    smy = smobj.y([smobj.ix(i) smobj.ixc(i)]);
    [xi,yi] = polyxpoly(d1(:,1),d1(:,2),smx,smy);
    if ~isempty(xi)
        irm = [irm; i];
    end
    [xi,yi] = polyxpoly(d2(:,1),d2(:,2),smx,smy);
    if ~isempty(xi)
        irm = [irm; i];
    end
end
smobj.ix(irm)  = [];
smobj.ixc(irm) = [];

bdi = [ [1:length(xbnd)-1]' [[2:length(xbnd)-1]';1] ones(length(xbnd)-1,1)];
smi = [length(bd) + smobj.ix length(bd) + smobj.ixc zeros(length(smobj.ix),1)];
% geom.point.coord = [bd;  sm ];
% geom.edge2.index = [bdi; smi];

geom.mshID = 'EUCLIDEAN-MESH';
if add_dam  
    i0  = length(bd);
    d1i = [ [i0+1:i0+length(d1)-1]' [i0+2:i0+length(d1)]' zeros(length(dam1.X)-2,1)];
    i0  = length(bd) + length(d1);
    d2i = [ [i0+1:i0+length(d2)-1]' [i0+2:i0+length(d2)]' zeros(length(dam2.X)-2,1)];
    i0  = length(bd) + length(d1) + length(d2);
    smi = [i0 + smobj.ix i0 + smobj.ixc zeros(length(smobj.ix),1)];
    geom.point.coord = [bd; d1;  d2 ];
    geom.edge2.index = [bdi;d1i; d2i];
else
    geom.point.coord = bd;
    geom.edge2.index = bdi;
end
if add_stream
    
end
figure(10);
coord = geom.point.coord;
index = geom.edge2.index;
for i = 1 : length(geom.edge2.index)
    plot([coord(index(i,1),1) coord(index(i,2),1)],[coord(index(i,1),2) coord(index(i,2),2)],'b-','LineWidth',2); hold on;
end

savemsh(opts.geom_file,geom) ;

opts.hfun_scal = 'absolute';
opts.hfun_hmax = +inf ;             % push HFUN limits
opts.mesh_dims = +2 ;               % 2-dim. simplexes
opts.optm_qlim = +.95 ;
opts.mesh_top1 = true ;             % for sharp feat's
%opts.mesh_top2 = true ; 
opts.geom_feat = true ;

hmat.mshID          = 'EUCLIDEAN-GRID';
hmat.point.coord{1} = xpos';
hmat.point.coord{2} = ypos ;
hmat.value          = hfun ;
hmat.slope          = dhdx*ones(size(hfun));
savemsh(opts.hfun_file,hmat);

hlim = marche(opts) ;

mesh = jigsaw(opts) ;

index_old = 1 : size(mesh.point.coord,1);
in = inpolygon(mesh.point.coord(:,1),mesh.point.coord(:,2),xbnd,ybnd);
mesh.point.coord(~in,:) = [];
index_old(~in) = [];
index_new = 1 : size(mesh.point.coord,1);

irm = [];
for i = 1 : size(mesh.tria3.index,1)
    if ~all(in(mesh.tria3.index(i,1:3)))
        irm = [irm; i];
    end
end
mesh.tria3.index(irm,:) = [];
% Reindex the edge2 and tria3
for i = 1 : length(index_old)
    for j = 1 : 3
        tmp = mesh.tria3.index(:,j);
        tmp(tmp == index_old(i)) = index_new(i);
        mesh.tria3.index(:,j) = tmp;
    end
    for j = 1 : 2
        tmp = mesh.edge2.index(:,j);
        tmp(tmp == index_old(i)) = index_new(i);
        mesh.edge2.index(:,j) = tmp;
    end
end
coordx  = mesh.point.coord(:,1);
coordy  = mesh.point.coord(:,2);
connect = mesh.tria3.index(:,1:3);
coordz  = griddata(xdem,ydem,dem,coordx,coordy);
if ~isempty(isnan(coordz))
    coordz(isnan(coordz)) = griddata(xdem,ydem,dem,coordx(isnan(coordz)),coordy(isnan(coordz)),'nearest');
end
coordz  = round(coordz,1);
coordb  = zeros(length(coordx),1);
[in,on] = inpolygon(coordx,coordy,xbnd,ybnd);
coordb(on) = 1;
ind_db  = find(coordb == 1);

%-#-#-#-#- IF DON NOT KNOW WHERE IS THE OUTLET -#-#-#-#-%
% ind_out = ind_db(find(coordz(ind_db) == min(coordz(ind_db))));
% circ = polybuffer([nanmean(coordx(ind_out)) nanmean(coordy(ind_out))],'points',2000);
% in = inpolygon(coordx(ind_db),coordy(ind_db),circ.Vertices(:,1),circ.Vertices(:,2));
% figure;
% plot(coordx(ind_db(in)),coordz(ind_db(in)),'ko','LineWidth',2); hold on;
% plot(coordx(ind_out),coordz(ind_out),'go','LineWidth',2); 
% circ = polybuffer([nanmean(coordx(ind_out)) nanmean(coordy(ind_out))],'points',rm_width);
%-#-#-#-#- *********************************** -#-#-#-#-%

fname2 = '../meshes/Turning_30m';
x2 = ncread([fname2 '.exo'],'coordx');
y2 = ncread([fname2 '.exo'],'coordy');
tri2 =  ncread([fname2 '.exo'],'connect1');
%triplot(double(tri2)',double(x2),double(y2)); hold on;
tmp = dlmread([fname2 '.bcode']);
bcode2 = tmp(2:end,4);
xlim([3.2315*1e6 3.235*1e6]);
ylim([3.633*1e6  3.637*1e6]);
plot(x2(bcode2 == 2),y2(bcode2 == 2),'gs','LineWidth',2);

circ = polybuffer([nanmean(x2(bcode2 == 2)) nanmean(y2(bcode2 == 2))],'points',200);
in = inpolygon(coordx(ind_db),coordy(ind_db),circ.Vertices(:,1),circ.Vertices(:,2));
ind_out = ind_db(in);
plot(coordx(ind_out),coordz(ind_out),'rx','LineWidth',2);

area = polyarea(coordx(connect'),coordy(connect'));
figure;
triplot(connect,coordx,coordy); hold on;
plot(dam1.X,dam1.Y,'g-','LineWidth',1);
plot(dam2.X,dam2.Y,'g-','LineWidth',1);
plot(coordx(ind_out),coordy(ind_out),'rx','LineWidth',2);
plot(circ);
coordb(ind_out) = 2;

jig2exo([rootpath name], mesh, coordb, coordz);


