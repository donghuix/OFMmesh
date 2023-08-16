function [gages,hwm,mesh] = read_outputs(gages,hwm,names,Output,ngrids)
    
    if nargin == 3
        ngrids = 200;
    end
    
    % Read mesh
    mesh = struct([]);
    [mesh(1).v_x,mesh(1).v_y,mesh(1).v_z,mesh(1).v_b] = plotbcode([],'../meshes/Turning_30m.bcode');
    mesh(1).tri = ncread('../meshes/Turning_30m.exo','connect1');

    nov = length(mesh(1).v_x)/2; 
    noc = size(mesh(1).tri,2); % number of cells
    mesh(1).c_x = mesh(1).v_x(mesh(1).tri);
    mesh(1).c_y = mesh(1).v_y(mesh(1).tri);
    mesh(1).c_b = mesh(1).v_b(mesh(1).tri);
    mesh(1).c_z = mesh(1).v_z(mesh(1).tri);

    mesh(1).x_tri = mean(mesh(1).c_x,1);
    mesh(1).y_tri = mean(mesh(1).c_y,1);
    mesh(1).z_tri = mean(mesh(1).c_z,1);

    % Identify index in mesh for the USGS and HWM 
    for i = 1 : length(gages)
        gages(i).idx = NaN(ngrids,1);
        dist  = (mesh.x_tri - gages(i).X).^2 + (mesh.y_tri - gages(i).Y).^2;
        [~,I] = sort(dist,'ascend');
        gages(i).idx = I(1:ngrids);
    end
    hwm.idx = NaN(length(hwm.X),ngrids);
    for i = 1 : length(hwm.X)
        dist  = (mesh.x_tri - hwm.X(i)).^2 + (mesh.y_tri - hwm.Y(i)).^2;
        [~,I] = sort(dist,'ascend');
        hwm.idx(i,:) = I(1:ngrids);
    end

    % Read OFM outputs
    t0 = datenum(2017,8,26,0,0,0);
    t1 = datenum(2017,8,30,23,59,59);
    tsim = t0 : 1/24 : t1; %tsim = tsim(1:72);
    tsim = tsim(1:end-1);

    nt = 119;
    for irun = 1 : 3
        name = names{irun};
        if ~exist([name '_Validations.mat'],'file')
            hgage = NaN(length(gages),ngrids,nt);
            hmark = NaN(length(hwm),ngrids,nt);
            for i = 1 : nt
                fnm = [Output name '/solution_' num2str(i) '.dat'];
                disp(fnm);
                arr = PetscBinaryRead(fnm);
                arr = reshape(arr, [3, noc]); % [h, uh, vh]
  
                for k = 1 : length(gages)
                    hgage(k,:,i) = arr(1,gages(k).idx);
                end
                for k = 1 : length(hwm.X)
                    hmark(k,:,i) = arr(1,hwm.idx(k,:));
                end
            end
            save([name '_Validations.mat'],'hgage','hmark');
        end
    end
    
    SR = load('Turning_30m_SR_Validations.mat');
    UF = load('Turning_30m_Validations.mat');
    NB = load('Turning_30m_noBC_SR_Validations.mat');
    delete('Turning_30m_Validations.mat');
    delete('Turning_30m_SR_Validations.mat');
    delete('Turning_30m_noBC_Validations.mat');

    % 
    hmark1 = max(SR.hmark,[],3); dist1 = abs(hmark1 - [hwm.h]);  
    hmark2 = max(UF.hmark,[],3); 
    hmark3 = max(NB.hmark,[],3); 
    for i = 1 : length(hwm.X)
        ind1 = find(dist1(i,:) == min(dist1(i,:)));
        hwm.sim(1).h(i,1) = hmark1(i,ind1);
        hwm.sim(2).h(i,1) = hmark2(i,ind1);
        hwm.sim(3).h(i,1) = hmark3(i,ind1);
    end

    % Read USGS observation
    for i = 1:length(gages)
        [twl,wl] = get_USGS_wl(gages(i).ID{1});
        wl = wl.*0.3048;
        if ~isempty(wl)
            gages(i).wl   = interp1(twl,wl,tsim)';
            gages(i).tsim = tsim;
            hsim1 = NaN(ngrids,nt);
            hsim2 = NaN(ngrids,nt);
            hsim3 = NaN(ngrids,nt);
            for j = 1 : ngrids
                hsim1(j,:) = SR.hgage(i,j,:)+min(mesh.c_z(:,gages(i).idx(j)));
                hsim2(j,:) = UF.hgage(i,j,:)+min(mesh.c_z(:,gages(i).idx(j)));
                hsim3(j,:) = NB.hgage(i,j,:)+min(mesh.c_z(:,gages(i).idx(j)));
                [R21(j),RMSE1(j),NSE1(j)] = estimate_evaluation_metric(gages(i).wl,hsim1(j,:)');
                [R22(j),RMSE2(j),NSE2(j)] = estimate_evaluation_metric(gages(i).wl,hsim2(j,:)');
                [R23(j),RMSE3(j),NSE3(j)] = estimate_evaluation_metric(gages(i).wl,hsim3(j,:)');
            end
            ind1 = find(NSE1 == max(NSE1));
            if length(ind1) > 1
                tmp = RMSE1(ind1);
                ind1 = ind1(tmp == min(tmp));
            end
            gages(i).sim(1).h    = hsim1(ind1,:);
            gages(i).sim(2).h    = hsim2(ind1,:);
            gages(i).sim(3).h    = hsim3(ind1,:);
            gages(i).sim(1).R2   = R21(ind1);
            gages(i).sim(2).R2   = R22(ind1);
            gages(i).sim(3).R2   = R23(ind1);
            gages(i).sim(1).RMSE = RMSE1(ind1);
            gages(i).sim(2).RMSE = RMSE2(ind1);
            gages(i).sim(3).RMSE = RMSE3(ind1);
            gages(i).sim(1).NSE  = NSE1(ind1);
            gages(i).sim(2).NSE  = NSE2(ind1);
            gages(i).sim(3).NSE  = NSE3(ind1);
            gages(i).sim(1).ex   = 'SR';
            gages(i).sim(2).ex   = 'UF';
            gages(i).sim(3).ex   = 'NB';
        else
            gages(i).wl          = [];
            gages(i).sim(1).h    = [];
            gages(i).sim(2).h    = [];
            gages(i).sim(3).h    = [];
            gages(i).sim(1).R2   = [];
            gages(i).sim(2).R2   = [];
            gages(i).sim(3).R2   = [];
            gages(i).sim(1).RMSE = [];
            gages(i).sim(2).RMSE = [];
            gages(i).sim(3).RMSE = [];
            gages(i).sim(1).NSE  = [];
            gages(i).sim(2).NSE  = [];
            gages(i).sim(3).NSE  = [];
            gages(i).sim(1).ex   = 'SR';
            gages(i).sim(2).ex   = 'UF';
            gages(i).sim(3).ex   = 'NB';
        end
    end

end