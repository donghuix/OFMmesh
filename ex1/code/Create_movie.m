clear;close all;clc;

[Output, exportfig] = SetupEnvironment();

load('../data/depth.mat','xbnd','ybnd','mesh');
proj = projcrs(32610);
[lat,lon] = projinv(proj,xbnd,ybnd);
[trilat,trilon] = projinv(proj,mesh.c_x,mesh.c_y);

cmap = getPanoply_cMap('NEO_omi_ozone_to3');
cmap(1:40,:) = [];

v1 = VideoWriter("Whole_domain.mp4",'MPEG-4');
v1.FrameRate = 5;
open(v1);

v2 = VideoWriter("Middle.mp4",'MPEG-4');
v2.FrameRate = 5;
open(v2);

v3 = VideoWriter("downscream.mp4",'MPEG-4');
v3.FrameRate = 5;
open(v3);

figure(1); set(gcf,'Position',[10 10 800 600]);
gx1 = geoaxes;
geoplot(gx1,lat,lon,'r-','LineWidth',2); hold on;
geoplot(gx1,[29.75 29.75 29.80 29.80 29.75],[-95.5 -95.4 -95.4 -95.5 -95.5],'b-','LineWidth',2);
geoplot(gx1,[29.75 29.75 29.775 29.775 29.75],[-95.375 -95.325 -95.325 -95.375 -95.375],'g-','LineWidth',2);
geobasemap(gx1,'satellite');
ax1 = axes;

figure(2); set(gcf,'Position',[810 10 800 600]);
gx2 = geoaxes;
geoplot(gx2,[29.75 29.75 29.80 29.80 29.75],[-95.5 -95.4 -95.4 -95.5 -95.5],'b-','LineWidth',2);
geobasemap(gx2,'satellite');
ax2 = axes;

figure(3); set(gcf,'Position',[810 610 800 600]);
gx3 = geoaxes;
geoplot(gx3,[29.75 29.75 29.775 29.775 29.75],[-95.375 -95.325 -95.325 -95.375 -95.375],'g-','LineWidth',2);
geobasemap(gx3,'satellite');
ax3 = axes;

t0 = datenum(2017,8,26,0,0,0);
for i = 1 : 119
    disp(i);
    fname = ['/Users/xudo627/Developments/ofm_petsc/Output/Turning_30m_SR/solution_' num2str(i) '.dat'];
    array = PetscBinaryRead(fname);
    array = reshape(array,[3 length(array)/3]);
    
    figure(1);
    delete(ax1);
    if i > 1
        delete(tit1);
    end
    ax1 = axes;
    ind = find(array(1,:) > 0.1);
    patch(trilon(:,ind),trilat(:,ind),array(1,ind),'LineStyle','none','FaceAlpha',0.7); colormap(cmap); clim([0 12]);
    xlim(gx1.LongitudeLimits);
    ylim(gx1.LatitudeLimits);
    ax1.Visible = 'off'; 
    ax1.XTick = []; 
    ax1.YTick = []; 
    t = t0 + i/24;
    [yr,mo,dy,hr] = datevec(t);
    str = [num2str(yr) '-0' num2str(mo) '-' num2str(dy) '  ' num2str(hr) ':00'];
    tit1 = add_title(ax1,str,30,'in'); tit1.Color = [1 1 1];
    
    frame1 = getframe(gcf);
    writeVideo(v1,frame1);

    figure(2);
    delete(ax2);
    if i > 1
        delete(tit2);
    end
    ax2 = axes;
    ind = find(array(1,:) > 0.1);
    patch(trilon(:,ind),trilat(:,ind),array(1,ind),'LineStyle','none','FaceAlpha',0.7); colormap(cmap); clim([0 12]);
    xlim(gx2.LongitudeLimits);
    ylim(gx2.LatitudeLimits);
    ax2.Visible = 'off'; 
    ax2.XTick = []; 
    ax2.YTick = []; 
    t = t0 + i/24;
    [yr,mo,dy,hr] = datevec(t);
    str = [num2str(yr) '-0' num2str(mo) '-' num2str(dy) '  ' num2str(hr) ':00'];
    tit2 = add_title(ax2,str,30,'in'); tit2.Color = [1 1 1];

    frame2 = getframe(gcf);
    writeVideo(v2,frame2);

    figure(3);
    delete(ax3);
    if i > 1
        delete(tit3);
    end
    ax3 = axes;
    ind = find(array(1,:) > 0.1);
    patch(trilon(:,ind),trilat(:,ind),array(1,ind),'LineStyle','none','FaceAlpha',0.7); colormap(cmap); clim([0 12]);
    xlim(gx3.LongitudeLimits - [0 0.001]);
    ylim(gx3.LatitudeLimits  - [0 0.002]);
    ax3.Visible = 'off'; 
    ax3.XTick = []; 
    ax3.YTick = []; 
    t = t0 + i/24;
    [yr,mo,dy,hr] = datevec(t);
    str = [num2str(yr) '-0' num2str(mo) '-' num2str(dy) '  ' num2str(hr) ':00'];
    tit3 = add_title(ax3,str,30,'in'); tit3.Color = [1 1 1];
    
    frame3 = getframe(gcf);
    writeVideo(v3,frame3);

    pause(0.1);
end
close(v1);
close(v2);
close(v3);