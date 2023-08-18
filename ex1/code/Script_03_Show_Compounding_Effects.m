clear;close all;clc;

[Output, exportfig] = SetupEnvironment();

labels = {'Less than Ankle','Ankle - Knee','Knee -Waist','Waist - Head','Above Head'};
% threshold for distance from a grid cell to outlet
thres  = 1000 : 1000 : 50000;

prefix = 'Turning_30m_SR';

if exist('../data/inundCF.mat','file')
    load('../data/inundCF.mat')
else
    [inundf,inundv,hdiff,inunds] = read_CF(Output, thres);
    save('../data/inundCF.mat','inunds','inundf','inundv','hdiff','-v7.3');
end

figure(4);

for i = 1 : 16
    subplot(4,4,i)
    tmp1 = sum(inundv(i,:,:,1),2);
    tmp2 = sum(inundv(i,:,:,2),2);
    plot(tmp1(:),'r-','LineWidth',2); hold on;
    plot(tmp2(:),'b-','LineWidth',2);
    plot(tmp1(:)-tmp2(:),'-','LineWidth',2); hold on; grid on;
end

figure(5);
tmp1 = sum(inundv(:,:,:,1),2);
tmp2 = sum(inundv(:,:,:,2),2);
tmp  = tmp2 - tmp1;





