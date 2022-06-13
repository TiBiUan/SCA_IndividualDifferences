%% X is the data at hand (n_samples x n_dims)
% idx contains the indexing information we will use for the GMM case (can
% be several)
function [Mu,Sigma,L_HC,LL_HC] = ET_FitGaussian(X)
    
    % Parameters
    Mu = mean(X)';
    Sigma = cov(X);

    % For control subjects, likelihood is computed for each region
    [L_HC,LL_HC] = ET_EvaluateGaussian(X',Mu,Sigma);

end