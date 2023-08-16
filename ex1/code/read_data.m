function [xbnd,ybnd,gages,hwm,dams] = read_data()

    % Projection coordinates
    proj = projcrs(32610);

    % Read boundary as longitude and latitude
    bndfile = '../turning.json';
    [xbnd,ybnd,dem,xdem,ydem] = getDEM(bndfile,'30m',false);

    % Read USGS gages
    fileID = fopen('../usgs/site_information');
    C = textscan(fileID,'%s%s%s%s%s%s%s%s%s%s%s%s','HeaderLines',33,'Delimiter','\t');
    fclose(fileID);
    gages = struct([]);
    k = 1;
    for i = 1 : length(C{9})
        if inpolygon(str2num(C{10}{i}),str2num(C{9}{i}),xbnd,ybnd)
            gages(k).Y  = str2num(C{9}{i}); 
            gages(k).X  = str2num(C{10}{i});
            gages(k).ID = C{1,2}(i);
            [gages(k).X,gages(k).Y] = projfwd(proj,gages(k).Y,gages(k).X);
            k = k + 1;
        end
    end

    % Read High Water Marks
    hwm = struct([]);
    T_hwm  = readtable('../hwm/Harvey_HWMs_20180322.csv');
    in_hwm = inpolygon(T_hwm.longitude,T_hwm.latitude,xbnd,ybnd);
    hwm(1).X = T_hwm.longitude(in_hwm);
    hwm(1).Y = T_hwm.latitude(in_hwm);
    hwm(1).h = T_hwm.height_above_gnd(in_hwm).*0.3048;
    hwm(1).X  = hwm(1).X(hwm(1).h > 0.001);
    hwm(1).Y  = hwm(1).Y(hwm(1).h > 0.001);
    hwm(1).h  = hwm(1).h(hwm(1).h > 0.001);
    [hwm(1).X,hwm(1).Y] = projfwd(proj,hwm(1).Y,hwm(1).X); % project to meter

    [xbnd,ybnd] = projfwd(proj,ybnd,xbnd);

    % add reservoir boundary
    tmp = shaperead('../dams/dam1.shp');
    dams(1).X = tmp.X + 0.00025;
    dams(1).Y = tmp.Y - 0.001;
    tmp = shaperead('../dams/dam2.shp');
    dams(2).X = tmp.X + 0.00036;
    dams(2).Y = tmp.Y - 0.00071;
    
    % Project dam coordinates
    [dams(1).X,dams(1).Y] = projfwd(proj,dams(1).Y,dams(1).X);
    [dams(2).X,dams(2).Y] = projfwd(proj,dams(2).Y,dams(2).X);
    ind = find(dams(2).Y < 3.627*1e6);
    dams(2).X(ind) = []; dams(2).Y(ind) = [];

end