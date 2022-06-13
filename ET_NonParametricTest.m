function [t,t_null,p,p_null] = ET_NonParametricTest(X1,X2,n_null)

    n1 = length(X1);
    n2 = length(X2);

    [~,p,~,stats] = ttest2(X1,X2);
    t = stats.tstat;
    
    for n = 1:n_null
        
        idx = randperm(n1+n2);
        tmp = [X1;X2];
        tmp = tmp(idx);
        
        X1n = tmp(1:n1);
        X2n = tmp(n1+1:end);
        
        [~,p_null(n),~,stats_null(n)] = ttest2(X1n,X2n);
        t_null(n) = stats_null(n).tstat;
    end
end