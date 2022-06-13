function [pmat,pvec,p_KLD] = ET_GetRegionProfile(RES,Alpha)

    n_perm = length(RES.DM1N);
    
    DM1T = prctile(RES.DM1N,[Alpha/2,100-Alpha/2]);
    DM2T = prctile(RES.DM2N,[Alpha/2,100-Alpha/2]);
    DM3T = prctile(RES.DM3N,[Alpha/2,100-Alpha/2]);
    
    DS11T = prctile(RES.DS11N,[Alpha/2,100-Alpha/2]);
    DS22T = prctile(RES.DS22N,[Alpha/2,100-Alpha/2]);
    DS33T = prctile(RES.DS33N,[Alpha/2,100-Alpha/2]);
    
    DS12T = prctile(RES.DS12N,[Alpha/2,100-Alpha/2]);
    DS13T = prctile(RES.DS13N,[Alpha/2,100-Alpha/2]);
    DS23T = prctile(RES.DS23N,[Alpha/2,100-Alpha/2]);
    
    p_KLD = min([2/n_perm*min([sum(RES.KLDN > RES.KLD),sum(RES.KLDN < RES.KLD)]),1]);
    
    p_DM1 = min([2/n_perm*min([sum(RES.DM1N > RES.DM1),sum(RES.DM1N < RES.DM1)]),1]);
    p_DM2 = min([2/n_perm*min([sum(RES.DM2N > RES.DM2),sum(RES.DM2N < RES.DM2)]),1]);
    p_DM3 = min([2/n_perm*min([sum(RES.DM3N > RES.DM3),sum(RES.DM3N < RES.DM3)]),1]);

    p_DS11 = min([2/n_perm*min([sum(RES.DS11N > RES.DS11),sum(RES.DS11N < RES.DS11)]),1]);
    p_DS22 = min([2/n_perm*min([sum(RES.DS22N > RES.DS22),sum(RES.DS22N < RES.DS22)]),1]);
    p_DS33 = min([2/n_perm*min([sum(RES.DS33N > RES.DS33),sum(RES.DS33N < RES.DS33)]),1]);

    p_DS12 = min([2/n_perm*min([sum(RES.DS12N > RES.DS12),sum(RES.DS12N < RES.DS12)]),1]);
    p_DS13 = min([2/n_perm*min([sum(RES.DS13N > RES.DS13),sum(RES.DS13N < RES.DS13)]),1]);
    p_DS23 = min([2/n_perm*min([sum(RES.DS23N > RES.DS23),sum(RES.DS23N < RES.DS23)]),1]);

    pmat = [p_DM1,p_DM2,p_DM3;p_DS11,p_DS22,p_DS33;p_DS12,p_DS13,p_DS23];
    
    pvec = [p_DM1;p_DM2;p_DM3;p_DS11;p_DS22;p_DS33;p_DS12;p_DS13;p_DS23];

end