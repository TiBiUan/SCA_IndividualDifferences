function [] = ET_PlotEllipses(Morpho_HC,Morpho_BASE,Morpho_YEAR,r,pro)

    X_HC = squeeze(Morpho_HC(r,:,:));
    [MU_HC,SIGMA_HC] = ET_FitGaussian(X_HC);

    X_BASE = squeeze(Morpho_BASE(r,:,:));
    [MU_BASE,SIGMA_BASE] = ET_FitGaussian(X_BASE);
    
    X_YEAR = squeeze(Morpho_YEAR(r,:,:));
    [MU_YEAR,SIGMA_YEAR] = ET_FitGaussian(X_YEAR);
    
    figure;
    
    plotcov3(MU_HC,SIGMA_HC,gca,[70,148,73]/255,0.1,pro);
    hold on
    plotcov3(MU_BASE,SIGMA_BASE,gca,[56,61,150]/255,0.1,pro);
    plotcov3(MU_YEAR,SIGMA_YEAR,gca,[175,54,60]/255,0.1,pro);
    scatter3(squeeze(Morpho_HC(r,:,1))',squeeze(Morpho_HC(r,:,2))',squeeze(Morpho_HC(r,:,3))',30,[157,188,64]/255,'filled');
    scatter3(squeeze(Morpho_BASE(r,:,1))',squeeze(Morpho_BASE(r,:,2))',squeeze(Morpho_BASE(r,:,3))',30,[8,133,161]/255,'filled');
    scatter3(squeeze(Morpho_YEAR(r,:,1))',squeeze(Morpho_YEAR(r,:,2))',squeeze(Morpho_YEAR(r,:,3))',30,[224,163,46]/255,'filled');

    xlabel('CT');
    ylabel('SA');
    zlabel('MC');






end