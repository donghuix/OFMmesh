clear;close all;clc;


load('boundary.mat');

debug = 0;

if strcmp(home,'/global/homes/d/donghui')
    %[nlcd,metadata] = readgeoraster('../data/nlcd_2021_land_cover_l48_20230630/nlcd_2021_land_cover_l48_20230630.img');
    
    lon = ncread('/global/cfs/projectdirs/m4267/donghui/grfr/RUNOFF_2017.nc','lon');
    lat = ncread('/global/cfs/projectdirs/m4267/donghui/grfr/RUNOFF_2017.nc','lat');
    
    disp(['lon dimension is ' num2str(length(lon))]);
    disp(['min(lon) = ' num2str(min(lon)) ', max(lon) = ' num2str(max(lon))]);
    disp(['lat dimension is ' num2str(length(lat))]);
    disp(['min(lat) = ' num2str(min(lat)) ', max(lat) = ' num2str(max(lat))]);
    
    %xmin = -125.75; xmax = -66.25; ymin = 24.25; ymax = 53.75;
    
    xmin = -115; xmax = -75; ymin = 27; ymax = 52;
    i1 = min(find(lon >= xmin));
    i2 = max(find(lon <= xmax));
    j1 = min(find(lat >= ymin));
    j2 = max(find(lat <= ymax));
    
    lon_region = ncread('/global/cfs/projectdirs/m4267/donghui/grfr/RUNOFF_2017.nc','lon',[i1],[i2-i1+1]);
    lat_region = ncread('/global/cfs/projectdirs/m4267/donghui/grfr/RUNOFF_2017.nc','lat',[j1],[j2-j1+1]);
    [lon_region,lat_region] = meshgrid(lon_region,lat_region);
    lon_region = lon_region'; lat_region = lat_region';
    
    in = inpoly2([lon_region(:) lat_region(:)],[xb_simp yb_simp]);
    ro = [];
    [m,n] = size(lon_region);
    
    for yr = 2017 : 2019
        filename = ['/global/cfs/projectdirs/m4267/donghui/grfr/RUNOFF_' num2str(yr) '.nc'];
        ro_region = NaN(m,n,365 * 8);
        for i = 1 : 365 * 8
            disp(['yr = ' num2str(yr) ', i = ' num2str(i)]);
            tmp = ncread(filename,'ro',[i1 j1 i],[i2-i1+1 j2-j1+1 1]);
            ro_region(:,:,i) = tmp;
            disp(['Runoff is ' num2str(nanmean(tmp(in)))]);
            ro = [ro; nanmean(tmp(in))];
        end
        save(['ro_region_' num2str(yr) '.mat'],'ro_region','lon_region','lat_region','-v7.3');
    end
    
    save('Mississippi_ro.mat','ro');
    
    figure;
    plot(ro,'k-','LineWidht',1); hold on; grid on;
    exportgraphics(gcf,'test.png','Resolution',300); 
    
    %imagesc(flipud(ro'));colorbar;colormap(jet);
    %caxis([0 5]);
    %exportgraphics(gcf,'test.png','Resolution',300); 
    
    
    if debug
        [lon,lat] = meshgrid(lon,lat);
        lon = lon';
        lat = lat';
        in = inpoly2([lon(:) lat(:)],[xb_simp yb_simp]);
        figure;
        plot(lon(in),lat(in),'rx'); hold on; grid on;
        plot(xb_simp,yb_simp,'k-','LineWidth',2); 
        exportgraphics(gcf,'test.png','Resolution',300);
    end

elseif strcmp(home,'/Users/xudo627')
    disp('Process Manning Coefficient!');
    mesh = load('mississippi.mat');
    
end


