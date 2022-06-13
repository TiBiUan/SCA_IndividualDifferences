function [R,R_null] = ET_NonParametricTest_corr(X1,X2,n_null)

    n1 = length(X1);
    n2 = length(X2);

    R = corr(X1,X2,'Type','Spearman');
    
    for n = 1:n_null
        
        idx = randperm(n1+n2);
        tmp = [X1;X2];
        tmp = tmp(idx);
        
        X1n = tmp(1:n1);
        X2n = tmp(n1+1:end);
        
        R_null(n) = corr(X1n,X2n,'Type','Spearman');
    end
end