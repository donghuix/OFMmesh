function [inundf,inundv,hdiff,inunds] = read_CF(Output, thres)

    fmesh  = '../meshes/Turning_30m';
    [x,y,z,b] = readbcode([fmesh '.bcode']);
    tri    = ncread([fmesh '.exo'],'connect1');
    % cell center
    xc = mean(x(tri)); yc = mean(y(tri)); zc = mean(z(tri),1);
    % cell vertices
    xv = x(tri);       yv = y(tri);
    % outlet
    xo = mean(x(b == 2)); yo = mean(y(b == 2)); 
    area = polyarea(xv, yv); noc  = length(zc);
    % distance from each cell to outlet
    dist = sqrt((xc-xo).^2 + (yc-yo).^2);
    clear tri;
    
    inunds = [0.01; 0.12; 0.46; 1; 1.70];
    % inund1 = find(arr(1,:) > 0.01 & arr(1,:) <= 0.12); % Less than Ankle
    % inund2 = find(arr(1,:) > 0.12 & arr(1,:) <= 0.46); % Ankle - Knee
    % inund3 = find(arr(1,:) > 0.46 & arr(1,:) <= 1.00); % Knee - Waist
    % inund4 = find(arr(1,:) > 1.00 & arr(1,:) <= 1.70); % Waist - Head
    % inund5 = find(arr(1,:) > 1.70 );                   % Above Head
    
    % Inundation Fraction
    inundf = NaN(length(thres),length(inunds),119,2);
    % Inundation Volume
    inundv = NaN(length(thres),119,2);
    
    for time = 1 : 119
        % Read simulated water depth [m] from outputs
        fnm1  = [Output 'Turning_30m_SR/solution_' num2str(time) '.dat'];
        arr1  = PetscBinaryRead(fnm1);
        arr1  = reshape(arr1, [3, noc]); % [h, uh, vh]
        hwtBC = arr1(1,:);
        %vmagwtBC = sqrt((arr1(2,:)./arr1(1,:)).^2 + (arr1(3,:)./arr1(1,:)).^2);
        
        fnm2  = [Output 'Turning_30m_noBC_SR/solution_' num2str(time) '.dat'];
        arr2  = PetscBinaryRead(fnm2);
        arr2  = reshape(arr2, [3, noc]); % [h, uh, vh]
        hnoBC = arr2(1,:);
        %vmagnoBC = sqrt((arr2(2,:)./arr2(1,:)).^2 + (arr2(3,:)./arr2(1,:)).^2);
        
        disp(fnm1);
        
        if time == 1
            hdiff = hwtBC' - hnoBC';
        else
            tmp   = hwtBC' - hnoBC';
            hdiff = max([hdiff tmp],[],2);
        end
    
        for i = 1 : length(thres)
            if i == 1
                ind = find(dist <= thres(i));
            else
                ind = find(dist <= thres(i) & dist > thres(i-1));
            end
    
            for j = 1 : length(inunds)
                tmph = hwtBC(ind);
                tmpa = area(ind);
                if j < 5
                    tmpi = find(tmph >= inunds(j) & tmph < inunds(j+1));
                else
                    tmpi = find(tmph > inunds(j));
                end
                
                inundf(i,j,time,1) = nansum(tmpa(tmpi))./nansum(tmpa);
                if j == 1
                    inundv(i,time,1) = nansum(tmph.*tmpa);
                end
        
                tmph = hnoBC(ind);
                tmpa = area(ind);
                if j < 5
                    tmpi = find(tmph >= inunds(j) & tmph < inunds(j+1));
                else
                    tmpi = find(tmph > inunds(j));
                end
                
                inundf(i,j,time,2) = nansum(tmpa(tmpi))./nansum(tmpa);
                if j == 1
                    inundv(i,time,2) = nansum(tmph.*tmpa);
                end
            end
        end
    end
    

end