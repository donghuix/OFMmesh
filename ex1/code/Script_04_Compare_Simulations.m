clear;close all;clc;

[Output, exportfig] = SetupEnvironment();

smeshes = {'Turning_30m', 'Turning_60m', 'Turning_90m', 'Turning_120m', ...
           'Turning_180m','Turning_240m','Turning_480m','Turning_990m'};

if exist('../data/ssims.mat','file')
    load('../data/ssims.mat');
else

    nt    = 119;
    ssims = struct([]);
    t0    = datenum(2017,8,26,0,0,0); % Start date
    
    for i = 1 : length(smeshes)
        mesh = [smeshes{i}];
        ssims(i).Name = mesh;
    
        % Read Discharge
        Qoutlet    = load([Output mesh '_SR/' mesh '_SR.Qoutlet']);
        ssims(i).t = Qoutlet(:,1)./86400 + t0;
        ssims(i).q = Qoutlet(:,2);
        
        % Read coordinates and connectivity from mesh file
        tri    = ncread(['../structure_meshes/' mesh '.exo'],'connect1');
        coordx = ncread(['../structure_meshes/' mesh '.exo'],'coordx');
        coordy = ncread(['../structure_meshes/' mesh '.exo'],'coordy');
        noc    = size(tri,2); % number of cells
        trix   = coordx(tri);
        triy   = coordy(tri); 
    
        % Read water depth and flood intensity and find the maximum values
        h  = zeros(noc,1);
        FI = zeros(noc,1);
        for j = 1 : nt
            fnm = [Output mesh '_SR/solution_' num2str(j) '.dat'];
            disp(['Reading ' mesh ': ' fnm]);
            arr = PetscBinaryRead(fnm);
            arr = reshape(arr, [3, noc]); % [h, uh, vh]
            tmp1 = arr(1,:);
            uh   = arr(2,:);
            vh   = arr(3,:);
            tmp2 = sqrt(uh.^2 + vh.^2);
    
            h    = max([h  tmp1'],[],2);
            FI   = max([FI tmp2'],[],2);
        end
        
        % interpolate to 30m meshes
        if i == 1
            trix0 = trix;
            triy0 = triy;
            load('MainChannel_poly.mat');
            in = inpoly2([mean(trix0)' mean(triy0)'],[polyout.Vertices(:,1) polyout.Vertices(:,2)]);
            ssims(i).hmax_wholedomain = h;
            ssims(i).hmax_mainchannel = ssims(i).hmax_wholedomain(in);
            ssims(i).FI_wholedomain   = FI;
            ssims(i).FI_mainchannel   = ssims(i).FI_wholedomain(in);
        else
            ssims(i).hmax_wholedomain = griddata(mean(trix,1)',mean(triy,1)',h, mean(trix0,1)',mean(triy0,1)');
            ssims(i).hmax_mainchannel = ssims(i).hmax_wholedomain(in);
            ssims(i).FI_wholedomain   = griddata(mean(trix,1)',mean(triy,1)',FI,mean(trix0,1)',mean(triy0,1)');
            ssims(i).FI_mainchannel   = ssims(i).FI_wholedomain(in);
        end
    end
    
    save('../data/ssims.mat','ssims');
end