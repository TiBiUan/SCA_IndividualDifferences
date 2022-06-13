function [Acc] = ET_SummonNaiveBayes(X1,X2,Mu1,Mu2,Sigma1,Sigma2)

    % Concatenates both populations and creates labels
    X = [X1;X2];
    y = [zeros(size(X1,1),1);ones(size(X2,1),1)];

    % Evaluates the result for each data point
    for s = 1:size(X,1)
        
        % Current data point
        x = X(s,:)';
        
        DeltaL(s) = ET_EvaluateGaussian(x,Mu1,Sigma1) - ET_EvaluateGaussian(x,Mu2,Sigma2);
    end
    
    y_pred = NaN(size(X,1),1);
    y_pred(DeltaL > 0) = 0;
    y_pred(DeltaL < 0) = 1;
    
    Acc = (size(X,1) - sum(abs(y - y_pred)))/size(X,1) * 100;
end