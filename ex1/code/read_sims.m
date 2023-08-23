function sims = read_sims(Output,meshes)
    
    nt    = 119;
    sims = struct([]);
    t0    = datenum(2017,8,26,0,0,0); % Start date
    
    fname = '../structure_meshes/Turning_30m.exo';
    tri0    = ncread(fname,'connect1');
    coordx = ncread(fname,'coordx');
    coordy = ncread(fname,'coordy');
    trix0   = coordx(tri0);
    triy0   = coordy(tri0); 
    clear coordx coordy tri0;
    load('MainChannel_poly.mat');
    in = inpoly2([mean(trix0)' mean(triy0)'],[polyout.Vertices(:,1) polyout.Vertices(:,2)]);

    for i = 1 : length(meshes)
        mesh = [meshes{i}];
        sims(i).Name = mesh;
    
        % Read Discharge
        Qoutlet    = load([Output mesh '_SR/' mesh '_SR.Qoutlet']);
        sims(i).t = Qoutlet(:,1)./86400 + t0;
        sims(i).q = Qoutlet(:,2);
        
        % Read coordinates and connectivity from mesh file
        fname = ['../meshes/' mesh '.exo'];
        if ~exist(fname,'file')
            fname = ['../structure_meshes/' mesh '.exo'];
        end
        if ~exist(fname,'file')
            error('Check the mesh dir!');
        end
        tri    = ncread(fname,'connect1');
        coordx = ncread(fname,'coordx');
        coordy = ncread(fname,'coordy');
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
        if size(trix,2) == size(trix0,2)
            sims(i).hmax_wholedomain = h;
            sims(i).hmax_mainchannel = sims(i).hmax_wholedomain(in);
            sims(i).FI_wholedomain   = FI;
            sims(i).FI_mainchannel   = sims(i).FI_wholedomain(in);
        else
            sims(i).hmax_wholedomain = griddata(mean(trix,1)',mean(triy,1)',h, mean(trix0,1)',mean(triy0,1)');
            sims(i).hmax_mainchannel = sims(i).hmax_wholedomain(in);
            sims(i).FI_wholedomain   = griddata(mean(trix,1)',mean(triy,1)',FI,mean(trix0,1)',mean(triy0,1)');
            sims(i).FI_mainchannel   = sims(i).FI_wholedomain(in);
        end
    end
        
end