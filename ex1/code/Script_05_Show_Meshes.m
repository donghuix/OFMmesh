clear;close all;clc;

[Output, exportfig] = SetupEnvironment();

meshs =  {'Turning_30_60_1500_1000_fpx2', ...
          'Turning_30_60_1500_1000',      ...
          'Turning_30_90_1500_1000',      ...
          'Turning_60_90_1500_1000',      ...
          'Turning_90_90_1500_1000',      ...
          'Turning_30_90_1500_10000',     ...
          'Turning_30_90_1500_100000',    ... 
          'Turning_30_90_1500_1000000',   ...             
          'Turning_30_90_1500_NoRiver',   ...
          'Turning_River_Dam_Burn_1150'};

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

figure;
for im = 1 : length(meshes)
    strs = strsplit(meshes{im},'_');
     % parameters to generate mesh
    hmain = str2double(strs{2});    % Main channle cell density
    htrib = str2double(strs{3});    % Tributary cell density 
    hhill = str2double(strs{4});    % Hillslope cell density
    raccu = str2double(strs{5});    % Accumulation threshold for river
    fp_width1  = 250;               % main channel floodplain width  [m]
    fp_width2  = 150;               % tributries floodplain width    [m]
    rm_width   = 180;               % river mouth width [m]  
    dhdx       = 1;                 % Marche gradient smooth
    add_dam    = 0;
    
    subplot(5,3,i);
    
    hfun  = ones(m,n).*hhill;
    
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
end

