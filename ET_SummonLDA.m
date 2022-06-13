function [Acc,Acc2] = ET_SummonLDA(X1,X2,Mu1,Mu2,Sigma1,Sigma2)

    % Concatenates both populations and creates labels
    X = [X1;X2];
    y = [zeros(size(X1,1),1);ones(size(X2,1),1)];

    Sigma = (Sigma1 + Sigma2)/2;
    w = Sigma * (Mu1 - Mu2);
    c = w' * 1/2 * (Mu1 + Mu2);
    
    % Evaluates the result for each data point
    for s = 1:size(X,1)
        
        % Current data point
        x = X(s,:)';
        
        DeltaL(s) = ((x - Mu1)' * inv(Sigma1) * (x - Mu1) + log(det(Sigma1))) - ...
            ((x - Mu2)' * inv(Sigma2) * (x - Mu2) + log(det(Sigma2)));
        
        DeltaL2(s) = w'*x - c;
    end
    
    y_pred = NaN(size(X,1),1);
    y_pred2 = y_pred;
    
    y_pred(DeltaL > 0) = 0;
    y_pred(DeltaL < 0) = 1;
    
    y_pred2(DeltaL2 > 0) = 0;
    y_pred2(DeltaL2 < 0) = 1;
    
    Acc = (size(X,1) - sum(abs(y - y_pred)))/size(X,1) * 100;
    Acc2 = (size(X,1) - sum(abs(y - y_pred2)))/size(X,1) * 100;
end