clear;close all;clc;

[Output, exportfig] = SetupEnvironment();

if exist('../data/rain.mat','file')
    load('../data/rain.mat');
else
    nt   = 119;
    rains = struct([]);
    t0   = datenum(2017,8,26,0,0,0); % Start date
    
    fname = '../structure_meshes/Turning_30m.exo';
    tri0    = ncread(fname,'connect1');
    coordx = ncread(fname,'coordx');
    coordy = ncread(fname,'coordy');
    trix0  = coordx(tri0);
    triy0  = coordy(tri0); 
    noc    = size(tri0,2);
    clear coordx coordy tri0;
    load('MainChannel_poly.mat');
    in = inpoly2([mean(trix0)' mean(triy0)'],[polyout.Vertices(:,1) polyout.Vertices(:,2)]);
    
    forcings = {'nldas','daymet'};
    
    for i = 1 : length(forcings)
        rains(i).Name = forcings{i};
    
        % Read Discharge
        Qoutlet    = load([Output 'Turning_30m_' forcings{i} '/' 'Turning_30m_' forcings{i} '.Qoutlet']);
        rains(i).t = Qoutlet(:,1)./86400 + t0;
        rains(i).q = Qoutlet(:,2);
    
        % Read water depth and flood intensity and find the maximum values
        h  = zeros(noc,1);
        FI = zeros(noc,1);
        for j = 1 : nt
            fnm = [Output 'Turning_30m_' forcings{i} '/solution_' num2str(j) '.dat'];
            disp(['Reading ' forcings{i} ': ' fnm]);
            arr = PetscBinaryRead(fnm);
            arr = reshape(arr, [3, noc]); % [h, uh, vh]
            tmp1 = arr(1,:);
            uh   = arr(2,:);
            vh   = arr(3,:);
            tmp2 = sqrt(uh.^2 + vh.^2);
    
            h    = max([h  tmp1'],[],2);
            FI   = max([FI tmp2'],[],2);
    
            rains(i).hmax_wholedomain = h;
            rains(i).hmax_mainchannel = rains(i).hmax_wholedomain(in);
            rains(i).FI_wholedomain   = FI;
            rains(i).FI_mainchannel   = rains(i).FI_wholedomain(in);
    
        end
    end
    
    save('../data/rain.mat','rains');
end