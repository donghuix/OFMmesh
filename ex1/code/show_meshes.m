clear;close all;clc;

addpath('../');
names  = {'Turning_30_60_1500_1000_fpx2','Turning_30_60_1500_1000', ...
          'Turning_30_90_1500_1000','Turning_60_90_1500_1000','Turning_90_90_1500_1000', ...
          'Turning_30_90_1500_10000','Turning_30_90_1500_100000','Turning_30_90_1500_1000000', ...              ...
          'Turning_30_90_1500_NoRiver','Turning_River_Dam_Burn_1150'};
names2  = {'Turning_30m','Turning_60m','Turning_90m','Turning_120m',...
           'Turning_180m','Turning_240m','Turning_480m','Turning_990m'};
figure(1);
for i = 1 : length(names2)
    name = names2{i};
    disp(name);
    if i == 1
        [x,y,z,b] = readbcode(['./meshes/' name '.bcode']);
        tri       = ncread(['./meshes/' name '.exo'],'connect1');
    else
        [x,y,z,b] = readbcode(['./structure_meshes/' name '.bcode']);
        tri       = ncread(['./structure_meshes/' name '.exo'],'connect1');
    end
    xtri = x(tri); ytri = y(tri); ztri = z(tri); 
    dz   = NaN(size(xtri,2),2);
    dz(:,1) = ((ytri(3,:) - ytri(1,:)) .* (ztri(2,:) - ztri(1,:)) -   ...
               (ytri(2,:) - ytri(1,:)) .* (ztri(3,:) - ztri(1,:))) ./ ...
              ((ytri(3,:) - ytri(1,:)) .* (xtri(2,:) - xtri(1,:)) -   ...
               (ytri(2,:) - ytri(1,:)) .* (xtri(3,:) - xtri(1,:)));

    dz(:,2) = ((xtri(3,:) - xtri(1,:)) .* (ztri(2,:) - ztri(1,:)) -   ...
               (xtri(2,:) - xtri(1,:)) .* (ztri(3,:) - ztri(1,:))) ./ ...
              ((ytri(3,:) - ytri(1,:)) .* (xtri(2,:) - xtri(1,:)) -   ...
               (ytri(2,:) - ytri(1,:)) .* (xtri(3,:) - xtri(1,:)));
    %[f,xi] = ecdf(sqrt(dz(:,1).^2 + dz(:,2).^2));
    [f,xi] = ksdensity(sqrt(dz(:,1).^2 + dz(:,2).^2),'Function','cdf');
    if i == 1
        semilogx(xi,f,'k-','LineWidth',2);hold on; grid on;
    else
        semilogx(xi,f,'-','LineWidth',2);
    end
end
legend(names2,'Interpreter','none');

for i = 1 : length(names)
    name = names{i};
    disp(name);
    [x,y,z,b] = readbcode(['./meshes/' name '.bcode']);
    tri       = ncread(['./meshes/' name '.exo'],'connect1');
    xtri = x(tri); ytri = y(tri); ztri = z(tri); 
    dz   = NaN(size(xtri,2),2);
    dz(:,1) = ((ytri(3,:) - ytri(1,:)) .* (ztri(2,:) - ztri(1,:)) -   ...
               (ytri(2,:) - ytri(1,:)) .* (ztri(3,:) - ztri(1,:))) ./ ...
              ((ytri(3,:) - ytri(1,:)) .* (xtri(2,:) - xtri(1,:)) -   ...
               (ytri(2,:) - ytri(1,:)) .* (xtri(3,:) - xtri(1,:)));

    dz(:,2) = ((xtri(3,:) - xtri(1,:)) .* (ztri(2,:) - ztri(1,:)) -   ...
               (xtri(2,:) - xtri(1,:)) .* (ztri(3,:) - ztri(1,:))) ./ ...
              ((ytri(3,:) - ytri(1,:)) .* (xtri(2,:) - xtri(1,:)) -   ...
               (ytri(2,:) - ytri(1,:)) .* (xtri(3,:) - xtri(1,:)));
    %[f,xi] = ecdf(sqrt(dz(:,1).^2 + dz(:,2).^2));
    [f,xi] = ksdensity(sqrt(dz(:,1).^2 + dz(:,2).^2),'Function','cdf');
    semilogx(xi,f,'--','LineWidth',2);
end
