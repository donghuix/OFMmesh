function plotmesh(coordx,coordy,coordz,coordb,connect)
    figure;
    patch(coordx(connect'),coordy(connect'),nanmean(coordz(connect),2),'LineStyle','none'); colorbar;
    hold on;
    plot(coordx(coordb == 2),coordy(coordb == 2),'kx','LineWidth',2);
end