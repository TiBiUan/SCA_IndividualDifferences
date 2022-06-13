% X contains the data to use for fitting (n_subjects x n_dims)
% K is the number of Gaussians to consider
function [MU,SIGMA,PI,LL] = ET_FitGMM_Classical(X,K,n_runs,epsilon)

    n_iter_max = 500;

    % Our dimensions
    n_dims = size(X,2);
    n_subjects = size(X,1);

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
            mu(:,k) = X(idx(k),:)';
            sigma(:,:,k) = cov(X(idx(1:n_subjects),:));
            ppi(k) = 1/K;
        end

        iter = 1;

        % To iterate until convergence
        while DL > epsilon && iter < n_iter_max
            
            % Expectation step
            for k = 1:K
                EVALS(:,k) = ppi(k)*ET_EvaluateGaussian(X',mu(:,k),squeeze(sigma(:,:,k)));
            end
            
            g = real(EVALS./(repmat(sum(EVALS,2),1,K)));

            % Maximization step
            ppi = real(mean(g)');

            mu = real((X' * g)./repmat(sum(g),n_dims,1));

            for k = 1:K
                sigma(:,:,k) = real(((repmat(g(:,k),1,n_dims).*(X - repmat(mu(:,k)',n_subjects,1)))'*(X - repmat(mu(:,k)',n_subjects,1)))./(realmin+sum(g(:,k))));
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
        clear g
        
        MU(:,:,n) = real(mu);
        SIGMA(:,:,:,n) = real(sigma);
        PI(:,n) = real(ppi);
        LL(n) = real(ll(end));
        
        clear sigma
        clear mu
        clear ppi
        clear ll
    end
end

