
function [twl,wl] = get_USGS_wl(site)

    verb = 0;

    fid = fopen('data');
    tline = fgetl(fid);
    k = 0;
    found_data = false;
    idx = 1;
    while ischar(tline)
        if contains(tline,['# Data provided for site ' site])
            disp(['Found ' site]);
            if verb == 1
                disp(tline);
            end
            k = 1;
        elseif contains(tline,'#            TS   parameter     Description') && k == 1
            tline = fgetl(fid);
            %disp(tline);
            if contains(tline,'00065')
                gageheight = true;
            end
        elseif contains(tline,['USGS	' site])
            found_data = true;
        elseif contains(tline,'#') && found_data
            k = 0;
        end
        if k > 0 && found_data
            if verb == 1
                disp(tline);
            end
            yr(idx) = str2num(tline(15:18));
            mo(idx) = str2num(tline(20:21));
            dy(idx) = str2num(tline(23:24));
            hr(idx) = str2num(tline(26:27));
            mi(idx) = str2num(tline(29:30));
            tmp     = strsplit(tline,'\t');
            twl(idx)= datenum(yr(idx),mo(idx),dy(idx),hr(idx),mi(idx),0);
            wl(idx) = str2num(tmp{5});
            x   = sscanf(tline,'%s %s %s %s %s %f %s');
            idx = idx + 1;
        end
        tline = fgetl(fid);
    end
    fclose(fid);

    if ~found_data
        twl = [];
        wl  = [];
    end
end