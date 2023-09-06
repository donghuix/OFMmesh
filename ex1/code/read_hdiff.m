clear;close all;clc;

[Output, exportfig] = SetupEnvironment();

nt   = 119;
sims = struct([]);
t0   = datenum(2017,8,26,0,0,0); % Start date

fname = '../structure_meshes/Turning_30m.exo';
tri0   = ncread(fname,'connect1');
coordx = ncread(fname,'coordx');
coordy = ncread(fname,'coordy');
trix0  = coordx(tri0);
triy0  = coordy(tri0); 
area   = polyarea(trix0, triy0)'; 
noc    = size(tri0,2);
clear coordx coordy tri0;
load('MainChannel_poly.mat');
in = inpoly2([mean(trix0)' mean(triy0)'],[polyout.Vertices(:,1) polyout.Vertices(:,2)]);

% Read water depth and flood intensity and find the maximum values
hdiff_BC      = zeros(noc,1);
hdiff_SR      = zeros(noc,1);
hdiff_BC_SR   = zeros(noc,1);

hdiff_BC2     = zeros(noc,1);
hdiff_SR2     = zeros(noc,1);
hdiff_BC_SR2  = zeros(noc,1);

hmax_BC_SR   = zeros(noc,1);
hmax_BC      = zeros(noc,1);
hmax_SR      = zeros(noc,1);
hmax         = zeros(noc,1);

inundv_BC_SR = zeros(nt,1);
inundv_BC    = zeros(nt,1);
inundv_SR    = zeros(nt,1);
inundv       = zeros(nt,1);


for j = 1 : nt
    fnm1 = [Output 'Turning_30m_SR/solution_' num2str(j) '.dat'];
    disp(['Reading: ' fnm1]);
    arr = PetscBinaryRead(fnm1);
    arr = reshape(arr, [3, noc]); % [h, uh, vh]
    h_BC_SR = arr(1,:);
    inundv_BC_SR(j) = nansum(h_BC_SR'.*area);
    
    fnm2 = [Output 'Turning_30m/solution_' num2str(j) '.dat'];
    disp(['Reading: ' fnm2]);
    arr = PetscBinaryRead(fnm2);
    arr = reshape(arr, [3, noc]); % [h, uh, vh]
    h_BC = arr(1,:);
    inundv_BC(j) = nansum(h_BC'.*area);

    fnm3 = [Output 'Turning_30m_noBC_SR/solution_' num2str(j) '.dat'];
    disp(['Reading: ' fnm3]);
    arr = PetscBinaryRead(fnm3);
    arr = reshape(arr, [3, noc]); % [h, uh, vh]
    h_SR = arr(1,:);
    inundv_SR(j) = nansum(h_SR'.*area);

    fnm4 = [Output 'Turning_30m_noBC/solution_' num2str(j) '.dat'];
    disp(['Reading: ' fnm4]);
    arr = PetscBinaryRead(fnm4);
    arr = reshape(arr, [3, noc]); % [h, uh, vh]
    h   = arr(1,:);
    inundv(j) = nansum(h'.*area);

    if j == 1
        hdiff_BC     = h_BC_SR' - h_SR';
        hdiff_SR     = h_BC_SR' - h_BC';
        hdiff_BC_SR  = h_BC_SR' - h';
        hdiff_BC2    = h_BC_SR' - h_SR';
        hdiff_SR2    = h_BC_SR' - h_BC';
        hdiff_BC_SR2 = h_BC_SR' - h';

        hmax_BC_SR  = h_BC_SR';
        hmax_BC     = h_BC';
        hmax_SR     = h_SR';
        hmax        = h';
    else

        hmax_BC_SR = max([hmax_BC_SR h_BC_SR'],[],2);
        hmax_BC    = max([hmax_BC    h_BC'   ],[],2);
        hmax_SR    = max([hmax_SR    h_SR'   ],[],2);
        hmax       = max([hmax       h'      ],[],2);

        tmp1 = h_BC_SR' - h_SR';
        tmp2 = h_BC_SR' - h_BC';
        tmp3 = h_BC_SR' - h';
        hdiff_BC    = max([hdiff_BC    tmp1],[],2);
        hdiff_SR    = max([hdiff_SR    tmp2],[],2);
        hdiff_BC_SR = max([hdiff_BC_SR tmp3],[],2);

        hdiff_BC2    = min([hdiff_BC2    tmp1],[],2);
        hdiff_SR2    = min([hdiff_SR2    tmp2],[],2);
        hdiff_BC_SR2 = min([hdiff_BC_SR2 tmp3],[],2);
%         b    = min([hdiff_BC    tmp1],[],2);
%         hdiff_BC    = a;
%         ind  = find(abs(b) > abs(a));
%         hdiff_BC(ind) = b(ind);

%         a    = max([hdiff_SR    tmp2],[],2);
%         b    = min([hdiff_SR    tmp2],[],2);
%         hdiff_SR    = a;
%         ind  = find(abs(b) > abs(a));
%         hdiff_SR(ind) = b(ind);
% 
%         a    = max([hdiff_BC_SR tmp3],[],2);
%         b    = min([hdiff_BC_SR tmp3],[],2);
%         hdiff_BC_SR = a;
%         ind  = find(abs(b) > abs(a));
%         hdiff_BC_SR(ind) = b(ind);

    end
end

save('../data/hdiff_decomp.mat','hdiff_BC','hdiff_SR','hdiff_BC_SR',     ...
                                'hdiff_BC2','hdiff_SR2','hdiff_BC_SR2',  ...
                                'hmax_BC_SR','hmax_BC','hmax_SR','hmax', ...
                                'inundv_BC_SR','inundv_BC','inundv_SR','inundv');