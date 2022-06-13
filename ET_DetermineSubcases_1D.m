%% This function determines which is the most likely data case at hand
function [Result] = ET_DetermineSubcases_1D(X,X_ref,K,n_runs,n_KLD,epsilon,alpha)

    n_subjects = length(X);

    % Number of parameters in both cases (Gaussian and GMM), used then for
    % the BIC computations
    D = 1;
    N = length(X);
    N_ref = length(X_ref);
    
    P_Gauss = 2*D + D*(D-1)/2;
    P_GMM = K*(2*D + D*(D-1)/2) + K;
    
    % Fits a Gaussian for the reference group and for the other (ML sense)
    [Mu_ref,Sigma_ref] = ET_FitGaussian_1D(X_ref);
    
    % Here, we want to loop and perform leave-one-out cross-validation
    for s = 1:n_subjects
        
        s
        
        % Segmenting into test and train data
        X_train = X;
        X_test = X(s);
        X_train(s) = [];
        
        % The model is fit on the training data
        [Mu_Gauss{s},Sigma_Gauss{s}] = ET_FitGaussian_1D(X_train');
        [Mu_GMM{s},Sigma_GMM{s},Pi_GMM{s},LL_GMM_train{s}] = ET_FitGMM_Classical_1D(X_train',K,n_runs,epsilon);
        idx_opt = find(LL_GMM_train{s} < Inf);
        idx_opt2 = find(LL_GMM_train{s}(idx_opt) == max(LL_GMM_train{s}(idx_opt)));
        Mu_GMM{s} = squeeze(Mu_GMM{s}(:,idx_opt(idx_opt2(1))));
        Sigma_GMM{s} = squeeze(Sigma_GMM{s}(:,idx_opt(idx_opt2(1))));
        Pi_GMM{s} = squeeze(Pi_GMM{s}(:,idx_opt(idx_opt2(1))));
        
        % The test data point is evaluated
        [~,LL_Gauss_test(s)] = ET_EvaluateGaussian_1D(X_test',Mu_Gauss{s},Sigma_Gauss{s});
        
        for k = 1:K
            EVALS(k) = real(Pi_GMM{s}(k)*ET_EvaluateGaussian_1D(X_test',Mu_GMM{s}(k),Sigma_GMM{s}(k)));
        end
        
        LL_GMM_test(s) = log(sum(EVALS));
    end

    % Total log-likelihood
    LL_Gauss = sum(LL_Gauss_test);
    LL_GMM = sum(LL_GMM_test);
    
    % Determines the Bayesian Information Criterion for both cases, and
    % also the Akaike Information Criterion
    AIC_Gauss = 2*P_Gauss - 2*LL_Gauss;
    AIC_GMM = 2*P_GMM - 2*LL_GMM;
    
    BIC_Gauss = P_Gauss*log(N) - 2*LL_Gauss;
    BIC_GMM = P_GMM*log(N) - 2*LL_GMM;
    
    % Decision time: Gaussian or GMM?
    if BIC_GMM < BIC_Gauss
        disp('GMM wins!');
    else
        disp('Gaussian wins!');
    end

    % We store all the parameters retrieved by the approaches
    Result.Mean_GMM = Mu_GMM;
    Result.Sigma_GMM = Sigma_GMM;
    Result.Pi_GMM = Pi_GMM;
    
    Result.Mean_Gauss = Mu_Gauss;
    Result.Sigma_Gauss = Sigma_Gauss;
    
    % We save the evaluation metrics as well
    Result.LL_GMM = LL_GMM;
    Result.LL_Gauss = LL_Gauss;
    
    Result.BIC_GMM = BIC_GMM;
    Result.BIC_Gauss = BIC_Gauss;
    
    Result.AIC_GMM = AIC_GMM;
    Result.AIC_Gauss = AIC_Gauss;
    
    [Mu_Gauss_all,Sigma_Gauss_all] = ET_FitGaussian_1D(X);
    
    % Kullback-Leibler divergence (actual)
    [KLD,DM,DS] = ET_ComputeExtendedGaussianStatistics_1D(Mu_ref,...
        Mu_Gauss_all,Sigma_ref,Sigma_Gauss_all);
        
    % We run many null settings to generate a null KLD distribution,
    % and compare it to the actual value. If the actual distance is
    % larger, then we know our second group differs from the reference
    for n = 1:n_KLD

        tmp_data = [X';X_ref'];
        id = randperm(N+N_ref);
        tmp_data = tmp_data(id);
        tmp_X = tmp_data(1:N);
        tmp_X_ref = tmp_data(N+1:end);

        [tmp_mu,tmp_sigma] = ET_FitGaussian_1D(tmp_X);
        [tmp_mu_ref,tmp_sigma_ref] = ET_FitGaussian_1D(tmp_X_ref);

        % Is there a significant group difference or not?
        [KLD_null(n),DM_null(n),DS_null(n)] = ...
            ET_ComputeExtendedGaussianStatistics_1D(tmp_mu_ref,tmp_mu,...
            tmp_sigma_ref,tmp_sigma);
    end
        
    KLD_thresh = prctile(KLD_null,100-alpha);
        
    % Saves the distance results
    Result.KLD = KLD;
    Result.KLDN = KLD_null;
    
    Result.DM = DM;
    Result.DMN = DM_null;
    
    Result.DS = DS;
    Result.DSN = DS_null;
        
    if KLD > KLD_thresh
        disp('Group difference (KLD)!');
        Result.GD = 1;
    else
        disp('No group difference!');
        Result.GD = 0;
    end

    % Which type of group difference do we have?
    % [Acc_QDA,Acc_LDA] = ET_SummonLDA(X,X_ref,Mu_Gauss_all,Mu_ref,Sigma_Gauss_all,Sigma_ref);
    Acc_QDA = NaN;
    Acc_LDA = NaN;
    Acc_NB = ET_SummonNaiveBayes_1D(X',X_ref',Mu_Gauss_all,Mu_ref,Sigma_Gauss_all,Sigma_ref);

    Result.Acc_QDA = Acc_QDA;
    Result.Acc_LDA = Acc_LDA;
    Result.Acc_NB = Acc_NB;
          
    Result.dis_LDA = NaN;
    Result.dis_NB = sqrt(sum((X' - repmat(Mu_Gauss_all,N,1)).^2,2));
end





