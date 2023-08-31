clear;close all;clc;

[Output, exportfig] = SetupEnvironment();

nt   = 119;
sims = struct([]);
t0   = datenum(2017,8,26,0,0,0); % Start date

fname = '../structure_meshes/Turning_30m.exo';
tri0    = ncread(fname,'connect1');
coordx = ncread(fname,'coordx');
coordy = ncread(fname,'coordy');
trix0   = coordx(tri0);
triy0   = coordy(tri0); 
noc    = size(tri0,2);
clear coordx coordy tri0;
load('MainChannel_poly.mat');
in = inpoly2([mean(trix0)' mean(triy0)'],[polyout.Vertices(:,1) polyout.Vertices(:,2)]);

% Read water depth and flood intensity and find the maximum values
hdiff_BC    = zeros(noc,1);
hdiff_SR    = zeros(noc,1);
hdiff_BC_SR = zeros(noc,1);

for j = 1 : nt
    fnm1 = [Output 'Turning_30m_SR/solution_' num2str(j) '.dat'];
    disp(['Reading: ' fnm1]);
    arr = PetscBinaryRead(fnm1);
    arr = reshape(arr, [3, noc]); % [h, uh, vh]
    h_BC_SR = arr(1,:);
    
    
    fnm2 = [Output 'Turning_30m/solution_' num2str(j) '.dat'];
    disp(['Reading: ' fnm2]);
    arr = PetscBinaryRead(fnm2);
    arr = reshape(arr, [3, noc]); % [h, uh, vh]
    h_BC = arr(1,:);

    fnm3 = [Output 'Turning_30m_noBC_SR/solution_' num2str(j) '.dat'];
    disp(['Reading: ' fnm3]);
    arr = PetscBinaryRead(fnm3);
    arr = reshape(arr, [3, noc]); % [h, uh, vh]
    h_SR = arr(1,:);

    fnm4 = [Output 'Turning_30m_noBC/solution_' num2str(j) '.dat'];
    disp(['Reading: ' fnm4]);
    arr = PetscBinaryRead(fnm4);
    arr = reshape(arr, [3, noc]); % [h, uh, vh]
    h   = arr(1,:);
    
    if j == 1
        hdiff_BC    = h_BC_SR' - h_SR';
        hdiff_SR    = h_BC_SR' - h_BC';
        hdiff_BC_SR = h_BC_SR' - h';
    else
        tmp1 = h_BC_SR' - h_SR';
        tmp2 = h_BC_SR' - h_BC';
        tmp3 = h_BC_SR' - h';
        hdiff_BC    = max([hdiff_BC    tmp1],[],2);
        hdiff_SR    = max([hdiff_SR    tmp2],[],2);
        hdiff_BC_SR = max([hdiff_BC_SR tmp3],[],2);
    end
end

save('../data/hdiff_decomp.mat','hdiff_BC','hdiff_SR','hdiff_BC_SR');