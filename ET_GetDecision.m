function [GD,LL_GMM,LL_Gauss,BIC_Gauss,BIC_GMM,AIC_Gauss,...
    AIC_GMM,Acc_NB,Acc_LDA,dis_NB,dis_LDA,KLD,KLDT,Pi_GMM] = ...
    ET_GetDecision(Result,n_regions,n_dims)
    
    switch n_dims
        
        case 1
            for i = 1:length(Result)
                GD(i) = Result{i}.GD;
                LL_GMM(i) = Result{i}.LL_GMM;
                LL_Gauss(i) = Result{i}.LL_Gauss;
                BIC_Gauss(i) = Result{i}.BIC_Gauss;
                BIC_GMM(i) = Result{i}.BIC_GMM;
                AIC_Gauss(i) = Result{i}.AIC_Gauss;
                AIC_GMM(i) = Result{i}.AIC_GMM;
                Acc_NB(i) = Result{i}.Acc_NB;
                Acc_LDA(i) = Result{i}.Acc_LDA;
                dis_NB(i,:) = Result{i}.dis_NB;
                dis_LDA(i,:) = Result{i}.dis_LDA;
                KLD(i) = Result{i}.KLD;
                KLDT(i) = prctile(Result{i}.KLDN,100-0.05/n_regions);
                Pi_GMM(i,:) = Result{i}.Pi_GMM;
            end
            
        case 2
            
            for i = 1:size(Result,1)
                for j = 1:size(Result,2)
                    GD(i,j) = Result{i,j}.GD;
                    LL_GMM(i,j) = Result{i,j}.LL_GMM;
                    LL_Gauss(i,j) = Result{i,j}.LL_Gauss;
                    BIC_Gauss(i,j) = Result{i,j}.BIC_Gauss;
                    BIC_GMM(i,j) = Result{i,j}.BIC_GMM;
                    AIC_Gauss(i,j) = Result{i,j}.AIC_Gauss;
                    AIC_GMM(i,j) = Result{i,j}.AIC_GMM;
                    Acc_NB(i,j) = Result{i,j}.Acc_NB;
                    Acc_LDA(i,j) = Result{i,j}.Acc_LDA;
                    dis_NB(i,j,:) = Result{i,j}.dis_NB;
                    dis_LDA(i,j,:) = Result{i,j}.dis_LDA;
                    KLD(i,j) = Result{i,j}.KLD;
                    KLDT(i,j) = prctile(Result{i,j}.KLDN,100-0.05/n_regions);
                    Pi_GMM(i,j,:) = Result{i,j}.Pi_GMM;
                end
            end
            
        case 3
            for i = 1:size(Result,1)
                for j = 1:size(Result,2)
                    for k = 1:size(Result,3)
                        GD(i,j,k) = Result{i,j,k}.GD;
                        LL_GMM(i,j,k) = Result{i,j,k}.LL_GMM;
                        LL_Gauss(i,j,k) = Result{i,j,k}.LL_Gauss;
                        BIC_Gauss(i,j,k) = Result{i,j,k}.BIC_Gauss;
                        BIC_GMM(i,j,k) = Result{i,j,k}.BIC_GMM;
                        AIC_Gauss(i,j,k) = Result{i,j,k}.AIC_Gauss;
                        AIC_GMM(i,j,k) = Result{i,j,k}.AIC_GMM;
                        Acc_NB(i,j,k) = Result{i,j,k}.Acc_NB;
                        Acc_LDA(i,j,k) = Result{i,j,k}.Acc_LDA;
                        dis_NB(i,j,k,:) = Result{i,j,k}.dis_NB;
                        dis_LDA(i,j,k,:) = Result{i,j,k}.dis_LDA;
                        KLD(i,j,k) = Result{i,j,k}.KLD;
                        KLDT(i,j,k) = prctile(Result{i,j,k}.KLDN,100-0.05/n_regions);
                        Pi_GMM(i,j,k,:) = Result{i,j,k}.Pi_GMM;
                    end
                end
            end
            
        otherwise
            errordlg('Error...');
    end
end