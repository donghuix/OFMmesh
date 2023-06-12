function [xbnd,ybnd,dem,xdem,ydem] = getDEM(bndfile,RES,showfigure)
% 
% This fucntion is to downlaod DEM for a given watershed boundary
%# Input:
% bndfile   : boundary file path + name, shape file, json file.
% RES       : resolution of DEM. '10m' or '30m'.
% showfigure: true or false.
[fundir,~,~] = fileparts(which('getDEM'));
if isempty(fundir)
    error('Cannot find getDEM! Need to compile OFMmesh first.');
end
[~, ~, fExt] = fileparts(bndfile);
if strcmp(fExt,'.shp')
    S = shaperead(bndfile);
    xbnd = S.X;
    ybnd = S.Y;
elseif strcmp(fExt,'.json')
    fid  = fopen(bndfile);
    raw  = fread(fid,inf);
    str  = char(raw');
    tmp  = jsondecode(str);
    xbnd = tmp.features.geometry.coordinates(1,:,1);
    ybnd = tmp.features.geometry.coordinates(1,:,2);
    fclose(fid);
end
xbnd = xbnd(:); ybnd = ybnd(:);
if ~exist('./DEM','dir')
    mkdir('./DEM');
end

fid   = fopen([fundir './DEM/NED_' RES '.txt']);
tline = fgetl(fid);
tiles = cell(2,1);
links = cell(2,1);
files = cell(2,1);
k     = 1;
xv    = NaN(4,2);
yv    = NaN(4,2);
while ischar(tline)
    strs     = strsplit(tline,'/');
    links{k} = tline;
    tiles{k} = strs{8};
    files{k} = strs{9};
    xv(1,k)  = -str2double(tiles{k}(5:7)) + 0;
    xv(2,k)  = -str2double(tiles{k}(5:7)) + 1;
    xv(3,k)  = -str2double(tiles{k}(5:7)) + 1;
    xv(4,k)  = -str2double(tiles{k}(5:7)) + 0;
    yv(1,k)  = +str2double(tiles{k}(2:3)) - 1;
    yv(2,k)  = +str2double(tiles{k}(2:3)) - 1;
    yv(3,k)  = +str2double(tiles{k}(2:3)) + 0;
    yv(4,k)  = +str2double(tiles{k}(2:3)) + 0;
    tline    = fgetl(fid);
    k        = k + 1;
end
fclose(fid);

ploybnd = polyshape(xbnd,ybnd);
idx     = false(length(tiles),1);
for i = 1 : length(tiles)
    polycell = polyshape(xv(:,i),yv(:,i));
    polyout  = intersect(ploybnd,polycell);
    if ~isempty(polyout.Vertices)
        idx(i) = true;
        if ~exist([fundir './DEM/' files{i}],'file')
            [status,cmdout] = system(['/usr/local/bin/wget ' links{i} ' -P ../DEM/'],'-echo');
            if status > 0
                disp(cmdout);
            end
        end
    end
end
idx = find(idx == 1);
idx = idx(2);

if showfigure
    figure;
    plot(xbnd,ybnd,'g-','LineWidth',4); hold on; grid on; axis equal;
    plot(xbnd,ybnd,'k-','LineWidth',1); 
    %plot(nanmean(xv,1),nanmean(yv,1),'k.');
end

for i = 1 : length(idx)
     polycell = polyshape(xv(:,idx(i)),yv(:,idx(i)));
     plot(polycell);
end

dem = double(imread([fundir '/../DEM/' files{idx}]));
dem(dem < -9999) = NaN;

[ny,nx] = size(dem);
dx = 1/nx;
dy = 1/ny;
xdem = xv(1,idx) + dx/2 : +dx : xv(2,idx) - dx/2;
ydem = yv(3,idx) - dy/2 : -dy : yv(1,idx) + dy/2;
[ydem,xdem] = meshgrid(ydem,xdem);

%in = inpolygon(xdem,ydem,xbnd,ybnd);

if showfigure
    figure;
    imagesc([xdem(1,1) xdem(end,end)],[ydem(1,1) ydem(end,end)],dem);
    set(gca,'YDir','normal'); colormap(gca,gray); hold on; 
    plot(xbnd,ybnd,'g-','LineWidth',4); grid on; axis equal;
    plot(xbnd,ybnd,'k-','LineWidth',1); 
end

end