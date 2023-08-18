clear;close all;clc;

[Output, exportfig] = SetupEnvironment();

labels = {'Less than Ankle','Ankle - Knee','Knee -Waist','Waist - Head','Above Head'};
% threshold for distance from a grid cell to outlet
thres  = 1000 : 1000 : 50000;

prefix = 'Turning_30m_SR';

if ~exist('../data/inundCF.mat','file')
    [inundf,inundv,hdiff,inunds] = read_CF(Output, thres);
    save('../data/inundCF.mat','inunds','inundf','inundv','hdiff','-v7.3');
end




