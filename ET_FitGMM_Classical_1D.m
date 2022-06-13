% X contains the data to use for fitting (n_subjects x n_dims)
% K is the number of Gaussians to consider
function [MU,SIGMA,PI,LL,GAMMA] = ET_FitGMM_Classical_1D(X,K,n_runs,epsilon)

    n_iter_max = 500;

    % Our dimensions
    n_subjects = length(X);

    % We will run many times with different parameter initializations
    for n = 1:n_runs

        % disp(['GMM run ',num2str(n),'...']);
        
        % This will be the difference in log-likelihood between successive
        % iterations
        DL = realmax;
        
        % Initialization, using random data points as starting values and
        % assuming the simplest possible covariance and pi values
        for k = 1:K
            idx = randperm(n_subjects);
            
            % Horizontal vectors
            mu(k) = X(idx(k));
            sigma(k) = std(X(idx(1:n_subjects)));
            ppi(k) = 1/K;
        end

        iter = 1;

        % To iterate until convergence
        while DL > epsilon && iter < n_iter_max
            
            % Expectation step
            for k = 1:K
                EVALS(:,k) = ppi(k)*ET_EvaluateGaussian_1D(X,mu(k),sigma(k));
            end
            
            g = real(EVALS./(repmat(sum(EVALS,2),1,K)));

            % Maximization step
            ppi = sum(g)/n_subjects;

            for k = 1:K
                mu(k) = (g(:,k)'*X)/sum(g(:,k));
                sigma(k) = sqrt(g(:,k)'*((X - mu(k)).^2)/sum(g(:,k)));
            end

            % Log-likelihood
            ll(iter) = sum(log(sum(EVALS,2)));

            if iter > 1
                DL = abs(ll(iter) - ll(iter-1));
            end

            if sum(isnan(sigma(:))) > 0
                break
            end
            
            iter = iter + 1;
        end
        
        clear EVALS
        
        MU(:,n) = real(mu);
        SIGMA(:,n) = real(sigma);
        PI(:,n) = real(ppi);
        LL(n) = real(ll(end));
        GAMMA(:,:,n) = g; 
        
        clear g
        clear sigma
        clear mu
        clear ppi
        clear ll
    end
end

