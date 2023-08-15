function [Output, proj] = SetupEnvironment()

    home = getenv('HOME');

    if strcmp(home,'/Users/xudo627')
        disp('Working on local!');

        addpath('/Users/xudo627/Developments/mylib/m/');
        addpath('/Users/xudo627/Developments/getPanoply_cMap/');
        % https://github.com/wschwanghart/topotoolbox
        addpath('/Users/xudo627/Developments/ofm_petsc/topotoolbox/');
        addpath('/Users/xudo627/Developments/ofm_petsc/topotoolbox/utilities/');
        addpath('/Users/xudo627/Developments/ofm_petsc/jigsaw-matlab/');
        addpath('/Users/xudo627/Developments/ofm_petsc/Matlab_Scripts/');
        % https://github.com/dengwirda/inpoly
        addpath('/Users/xudo627/Developments/inpoly/');
        % https://github.com/donghuix/Setup-E3SM-Mac
        addpath('/Users/xudo627/Developments/Setup-E3SM-Mac/matlab-scripts-to-process-inputs/');
        % https://github.com/g2e/m_map
        addpath('/Users/xudo627/Developments/m_map/');
        Output = '/Users/xudo627/Developments/ofm_petsc/Output/';

    elseif strcmp(home,'/qfs/people/xudo627')
        disp('Working on compy!');
        
        addpath('/qfs/people/xudo627/mylib/m/');
        addpath('/qfs/people/xudo627/getPanoply_cMap/');
        addpath('/compyfs/xudo627/ofm_petsc/topotoolbox/');
        addpath('/compyfs/xudo627/ofm_petsc/topotoolbox/utilities/');
        addpath('/compyfs/xudo627/ofm_petsc/jigsaw-matlab/');
        addpath('/compyfs/xudo627/ofm_petsc/Matlab_Scripts/');
        addpath('/qfs/people/xudo627/Setup-E3SM-Mac/matlab-scripts-to-process-inputs/');
        addpath('/qfs/people/xudo627/inpoly/');
        addpath('/qfs/people/xudo627/m_map/');
        Output = '/compyfs/xudo627/ofm_petsc/Output/';

    else
        disp('Need to install necessary packages!')
        stop(['Cannot work on ' home '!']);

    end

    addpath('../usgs/');
    addpath('../../');

    % Projection coordinates
    proj = projcrs(32610);
end