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

load('../data/ssims.mat');

load('colorblind_colormap.mat');
colorblind(3,:) = colorblind(12,:);
colorblind(9,:) = colorblind(11,:);

figure;
plot(ssims(1).t,ssims(1).q,'k-','LineWidth',2); hold on; grid on;
for i = 1 : length(rains)
    plot(rains(i).t,rains(i).q,'-','Color',colorblind(i+1,:),'LineWidth',2); 
end
datetick('x','keeplimits');
ylabel('Discharge [m^{3}/s]','FontSize',15,'FontWeight','bold');
legend('MRMS','NLDAS','daymet');

tri    = ncread('../structure_meshes/Turning_30m.exo','connect1');
coordx = ncread('../structure_meshes/Turning_30m.exo','coordx');
coordy = ncread('../structure_meshes/Turning_30m.exo','coordy');
trix   = coordx(tri);
triy   = coordy(tri); 
cmap = getPanoply_cMap('NEO_modis_lst');

figure;
ax(1+i) = subplot(1,3,1);
patch(trix,triy,ssims(1).hmax_wholedomain,'LineStyle','none'); clim([0 12]);
colormap(gca,cmap);
for i = 1 : 2
    ax(1+i) = subplot(1,3,i+1);
    patch(trix,triy,rains(i).hmax_wholedomain,'LineStyle','none'); clim([0 12]);
    set(gca,'XTick',[],'YTick',[]); colormap(gca,cmap);

    [R2(i),RMSE(i),~,PBIAS(i)] = estimate_evaluation_metric(ssims(1).hmax_wholedomain,rains(i).hmax_wholedomain);

    if i == 8
        cb = colorbar('south');
    end
end



figure;
histogram(ssims(1).hmax_wholedomain); hold on;
histogram(rains(1).hmax_wholedomain);