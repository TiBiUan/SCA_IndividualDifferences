function [pvec,p_KLD] = ET_GetRegionProfile_1D(RES)

    n_perm = length(RES.DMN);

    p_KLD = min([2/n_perm*min([sum(RES.KLDN > RES.KLD),sum(RES.KLDN < RES.KLD)]),1]);
    
    p_DM = min([2/n_perm*min([sum(RES.DMN > RES.DM),sum(RES.DMN < RES.DM)]),1]);
    
    p_DS = min([2/n_perm*min([sum(RES.DSN > RES.DS),sum(RES.DSN < RES.DS)]),1]);

    pvec = [p_DM;p_DS];

end