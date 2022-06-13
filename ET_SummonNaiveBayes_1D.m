function [Acc] = ET_SummonNaiveBayes_1D(X1,X2,Mu1,Mu2,Sigma1,Sigma2)

    % Concatenates both populations and creates labels
    X = [X1;X2];
    y = [zeros(length(X1),1);ones(length(X2),1)];

    % Evaluates the result for each data point
    for s = 1:length(X)
        
        % Current data point
        x = X(s);
        
        DeltaL(s) = ET_EvaluateGaussian_1D(x,Mu1,Sigma1) - ET_EvaluateGaussian_1D(x,Mu2,Sigma2);
    end
    
    y_pred = NaN(length(X),1);
    y_pred(DeltaL > 0) = 0;
    y_pred(DeltaL < 0) = 1;
    
    Acc = (length(X) - sum(abs(y - y_pred)))/length(X) * 100;
    
end