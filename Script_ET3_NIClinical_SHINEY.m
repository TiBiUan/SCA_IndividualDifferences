%% This script introduces an approach to analyse morphometric data while
% accounting for the dependences between different properties; i.e., it
% offers a tool for multivariate morphometric analyses. The framework
% relies on simple Bayesian probabilistic concepts.
% The approach is applied to a set of patients with essential
% tremor, who are compared to matched healthy controls.
% The framework reveals group differences in variance, and two different 
% sets of association between morphometric data in ET patients and clinical 
% symptoms.



%% Part 1 Loading of the data and generation of colormaps

addpath(genpath('Data'));
addpath(genpath('Utilities'));

% Load the data from all 3 groups
load ET_Clinical
load ET_Freesurfer
load Age_HC
load Gender_HC
load ET_Covariates
load CodeBook
load('FamilyHistory.mat');
load ET_Freesurfer_SC

% Decides the type of metric to probe
Metric = 'Thickness';
Metric2 = 'Area';
Metric3 = 'MeanCurvature';

% Cortical thickness
X_HC = (eval(['HC_',Metric]))';
X_BASE = (eval(['BASE_',Metric]))';
X_YEAR = (eval(['YEAR_',Metric]))';

% Surface area
Y_HC = (eval(['HC_',Metric2]))';
Y_BASE = (eval(['BASE_',Metric2]))';
Y_YEAR = (eval(['YEAR_',Metric2]))';

% Surface area
Z_HC = (eval(['HC_',Metric3]))';
Z_BASE = (eval(['BASE_',Metric3]))';
Z_YEAR = (eval(['YEAR_',Metric3]))';

% Number of subjects in total
n_HC = size(X_HC,2);
n_BASE = size(X_BASE,2);
n_YEAR = size(X_YEAR,2);

% Number of regions in the parcellation at play
n_regions = size(X_HC,1);
n_regions_SC = size(HC_Volume_SC,2);

% Will we plot the data or not in the script?
is_plot = 1;

% Tolerance level for convergence of log-likelihood
tol = 1e-5;


%%%%%%%% COLORMAPS

% Colormap to use for some representations (red-blue)
CM_RB = cbrewer('div','RdBu',1001);
CM_RB(CM_RB<0) = 0;

CM_RYG = cbrewer('div','RdYlGn',1001);
CM_RYG(CM_RYG < 0) = 0;

CM_Paired = cbrewer('qual','Paired',12);
CM_Paired(CM_Paired < 0) = 0;

% Mapping to lobes
LobesMapping = [[5,9,1,7,5,5,3,5,9,7,1,7,1,5,5,1,1,1,1,7,3,9,1,3,9,1,...
    1,3,5,3,1,5,5,11],1+[5,9,1,7,5,5,3,5,9,7,1,7,1,5,5,1,1,1,1,7,3,9,1,3,9,1,1,3,5,3,1,5,5,11]]';

CM_YG = cbrewer('seq','YlGn',1001);
CM_YG(CM_YG < 0) = 0;

CM_Subjects = colorcube(n_BASE);
CM_Subjects2 = flipud(cbrewer('div','RdYlGn',n_BASE));
CM_Subjects2(CM_Subjects2 < 0) = 0;

CM_tmp = cbrewer('qual','Set2',5);
CM_tmp1 = cbrewer('seq','Reds',1000);
CM_tmp1 = CM_tmp1(200:end,:);

CM_tmp2 = cbrewer('seq','Greens',1000);
CM_tmp2 = CM_tmp2(200:end,:);



%% Part 2 Regressing out covariates from the data
% Here, the data is regressed out for confounding variables: age, gender
% and total grey matter volume. In order to account for the fact that we
% considered two related groups (two measurements per subject, one in the
% pre-intervention and one in the post-intervention groups), we use a mixed 
% model approach instead of a group-wise regression. 
% Notwithstanding the within-subject variance modelling, this approach is 
% identical to group-wise regression if the three confounding factors are 
% considered jointly with their interaction with the group factor. 
% We regress out covariates of no interest, and keep the residuals as well 
% as the within-subject effects

% Data for mixed modelling
Age_MM = [Age_HC;Age;Age+1];
Gender_MM = [Gender_HC;Gender;Gender];
Subject_MM = nominal([((n_BASE+1:n_BASE+n_HC)');((1:n_BASE)');((1:n_YEAR)')]);
TGV_MM = [HC_TGV;BASE_TGV;YEAR_TGV];
Group_MM = nominal([ones(n_HC,1);2*ones(n_BASE,1);3*ones(n_YEAR,1)]);

% We perform the same process for each region...
for r = 1:n_regions
    
    % Morphometric data at hand...
    Morpho_MM = [X_HC(r,:)';X_BASE(r,:)';X_YEAR(r,:)'];

    % All summarized in a table to enable mixed modelling...
    TBL = table(Morpho_MM,Age_MM,Gender_MM,TGV_MM,Subject_MM,Group_MM);

    % Mixed model itself
    lme = fitlme(TBL,'Morpho_MM ~ 1 + Group_MM*(Age_MM + TGV_MM + Gender_MM) - Group_MM + (1 | Subject_MM)',...
        'FitMethod','ML');

    % The residuals are considered
    RES = lme.residuals;
    
    X_HC_res(r,:) = RES(1:n_HC);
    X_BASE_res(r,:) = RES(n_HC+1:n_HC+n_BASE);
    X_YEAR_res(r,:) = RES(n_HC+n_BASE+1:end);
    
    X_HC_res(r,:) = X_HC_res(r,:) + lme.Coefficients.Estimate(1);
    X_BASE_res(r,:) = X_BASE_res(r,:) + lme.Coefficients.Estimate(1);
    X_YEAR_res(r,:) = X_YEAR_res(r,:) + lme.Coefficients.Estimate(1);
    
    clear TBL
    
    % Same for the other morphometric properties
    Morpho_MM = [Y_HC(r,:)';Y_BASE(r,:)';Y_YEAR(r,:)'];

    TBL = table(Morpho_MM,Age_MM,Gender_MM,TGV_MM,Subject_MM,Group_MM);

    lme = fitlme(TBL,'Morpho_MM ~ 1 + Group_MM*(Age_MM + TGV_MM + Gender_MM) - Group_MM + (1 | Subject_MM)',...
        'FitMethod','ML');

    RES = lme.residuals;
    
    Y_HC_res(r,:) = RES(1:n_HC);
    Y_BASE_res(r,:) = RES(n_HC+1:n_HC+n_BASE);
    Y_YEAR_res(r,:) = RES(n_HC+n_BASE+1:end);
    
    Y_HC_res(r,:) = Y_HC_res(r,:) + lme.Coefficients.Estimate(1);
    Y_BASE_res(r,:) = Y_BASE_res(r,:) + lme.Coefficients.Estimate(1);
    Y_YEAR_res(r,:) = Y_YEAR_res(r,:) + lme.Coefficients.Estimate(1);
    
    clear TBL
    
    Morpho_MM = [Z_HC(r,:)';Z_BASE(r,:)';Z_YEAR(r,:)'];

    TBL = table(Morpho_MM,Age_MM,Gender_MM,TGV_MM,Subject_MM,Group_MM);

    lme = fitlme(TBL,'Morpho_MM ~ 1 + Group_MM*(Age_MM + TGV_MM + Gender_MM) - Group_MM + (1 | Subject_MM)',...
        'FitMethod','ML');

    RES = lme.residuals;
    
    Z_HC_res(r,:) = RES(1:n_HC);
    Z_BASE_res(r,:) = RES(n_HC+1:n_HC+n_BASE);
    Z_YEAR_res(r,:) = RES(n_HC+n_BASE+1:end);
    
    Z_HC_res(r,:) = Z_HC_res(r,:) + lme.Coefficients.Estimate(1);
    Z_BASE_res(r,:) = Z_BASE_res(r,:) + lme.Coefficients.Estimate(1);
    Z_YEAR_res(r,:) = Z_YEAR_res(r,:) + lme.Coefficients.Estimate(1);
    
    clear TBL
    
    if r <= n_regions_SC
        Morpho_MM = [HC_Volume_SC(:,r);BASE_Volume_SC(:,r);YEAR_Volume_SC(:,r)];

        TBL = table(Morpho_MM,Age_MM,Gender_MM,TGV_MM,Subject_MM,Group_MM);

        lme = fitlme(TBL,'Morpho_MM ~ 1 + Group_MM*(Age_MM + TGV_MM + Gender_MM) - Group_MM + (1 | Subject_MM)',...
            'FitMethod','ML');

        RES = lme.residuals;
        
        SC_HC_res(r,:) = RES(1:n_HC);
        SC_BASE_res(r,:) = RES(n_HC+1:n_HC+n_BASE);
        SC_YEAR_res(r,:) = RES(n_HC+n_BASE+1:end);
        
        SC_HC_res(r,:) = SC_HC_res(r,:) + lme.Coefficients.Estimate(1);
        SC_BASE_res(r,:) = SC_BASE_res(r,:) + lme.Coefficients.Estimate(1);
        SC_YEAR_res(r,:) = SC_YEAR_res(r,:) + lme.Coefficients.Estimate(1);
        
        clear TBL
    end
end

% Organizing the data more properly
Morpho_HC(:,:,1) = X_HC_res;
Morpho_HC(:,:,2) = Y_HC_res;
Morpho_HC(:,:,3) = Z_HC_res;

Morpho_BASE(:,:,1) = X_BASE_res;
Morpho_BASE(:,:,2) = Y_BASE_res;
Morpho_BASE(:,:,3) = Z_BASE_res;

Morpho_YEAR(:,:,1) = X_YEAR_res;
Morpho_YEAR(:,:,2) = Y_YEAR_res;
Morpho_YEAR(:,:,3) = Z_YEAR_res;

Morpho_BASE_SC = SC_BASE_res;
Morpho_YEAR_SC = SC_YEAR_res;
Morpho_HC_SC = SC_HC_res;



%% Part 3. Analysis of actual data (HC vs ET pre-intervention, and pre-
% versus post-intervention) in a "classical" fashion

% Symptoms before the intervention: we use everything we have to our
% disposal
Symptoms_pre = [ADL_Baseline,TSTH_Baseline,HEAD_Baseline,SymptomsDuration,FamilyHistory];

%%%%% "Classical" approach: we assume homogeneous groups 
% (non-clustered data). This step essentially treats the three properties
% independently from each other, hence likely missing on some important
% information. We correct for the number of tests performed in parallel
% Note that we perform t-testing; thus, we assume that the data is normally
% distributed (although we assess significance non-parametrically via
% permutation testing)

% Cortical data (3 properties)
for r = 1:n_regions
    
    r
    
    [t_CT_HB(r),t_null_CT_HB(r,:)] = ET_NonParametricTest(squeeze(Morpho_HC(r,:,1))',squeeze(Morpho_BASE(r,:,1))',10000);
    [t_SA_HB(r),t_null_SA_HB(r,:)] = ET_NonParametricTest(squeeze(Morpho_HC(r,:,2))',squeeze(Morpho_BASE(r,:,2))',10000);
    [t_MC_HB(r),t_null_MC_HB(r,:)] = ET_NonParametricTest(squeeze(Morpho_HC(r,:,3))',squeeze(Morpho_BASE(r,:,3))',10000);
end

% Non-cortical data (only volume)
for r = 1:n_regions_SC
    
    r
    
    [t_VOL_HB(r),t_null_VOL_HB(r,:)] = ET_NonParametricTest(Morpho_HC_SC(r,:)',Morpho_BASE_SC(r,:)',10000);
end

% Bonferroni-corrected significance threshold
t_CT_HB_thresh = prctile(t_null_CT_HB,[2.5/(3*n_regions + n_regions_SC),100-2.5/(3*n_regions + n_regions_SC)],2);
t_SA_HB_thresh = prctile(t_null_SA_HB,[2.5/(3*n_regions + n_regions_SC),100-2.5/(3*n_regions + n_regions_SC)],2);
t_MC_HB_thresh = prctile(t_null_MC_HB,[2.5/(3*n_regions + n_regions_SC),100-2.5/(3*n_regions + n_regions_SC)],2);

t_VOL_HB_thresh = prctile(t_null_VOL_HB,[2.5/(3*n_regions + n_regions_SC),100-2.5/(3*n_regions + n_regions_SC)],2);
t_VOL_YB_thresh = prctile(t_null_VOL_YB,[2.5/(3*n_regions + n_regions_SC),100-2.5/(3*n_regions + n_regions_SC)],2);

idx_CT_HB = find(t_CT_HB' < t_CT_HB_thresh(:,1) | t_CT_HB' > t_CT_HB_thresh(:,2));
idx_SA_HB = find(t_SA_HB' < t_SA_HB_thresh(:,1) | t_SA_HB' > t_SA_HB_thresh(:,2));
idx_MC_HB = find(t_MC_HB' < t_MC_HB_thresh(:,1) | t_MC_HB' > t_MC_HB_thresh(:,2));
idx_VOL_HB = find(t_VOL_HB' < t_VOL_HB_thresh(:,1) | t_VOL_HB' > t_VOL_HB_thresh(:,2));

% There is essentially nothing significant out of these above tests... This
% means that our data does not showcase any difference in terms of mean
% values!

% Association between pre-intervention data and symptoms: link between MC
% in largest part and HEAD + Family History
X0 = [Morpho_BASE(:,:,1)',Morpho_BASE(:,:,2)',Morpho_BASE(:,:,3)',Morpho_BASE_SC'];
Y0 = Symptoms_pre;
label = ones(n_BASE,1);
myPLS_inputs;
[IN_classical_pre,pls_opts_classical_pre] = myPLS_initialize(input,pls_opts,save_opts);
PLS_classical_pre = myPLS_analysis(IN_classical_pre,pls_opts_classical_pre);

if is_plot
    
    figure;
    bar(PLS_classical_pre.boot_results.Vb_mean(:,1));
    hold on
    errorbar(1:3*n_regions+n_regions_SC,PLS_classical_pre.boot_results.Vb_mean(:,1),...
        PLS_classical_pre.boot_results.Vb_mean(:,1)-PLS_classical_pre.boot_results.Vb_lB(:,1),...
        PLS_classical_pre.boot_results.Vb_uB(:,1)-PLS_classical_pre.boot_results.Vb_mean(:,1),'color','k','LineStyle','None','CapSize',3);

    figure;
    bar(PLS_classical_pre.boot_results.Ub_mean(:,1));
    hold on
    errorbar(1:size(Symptoms_pre,2),PLS_classical_pre.boot_results.Ub_mean(:,1),...
        PLS_classical_pre.boot_results.Ub_mean(:,1)-PLS_classical_pre.boot_results.Ub_lB(:,1),...
        PLS_classical_pre.boot_results.Ub_uB(:,1)-PLS_classical_pre.boot_results.Ub_mean(:,1),'color','k','LineStyle','None','CapSize',3);
end



%% Part 4. Comparison of Gaussian vs GMM model to characterize the areas
% and assessment of group differences

% We assume two clusters for now (i.e., consider the possibility of bimodal
% distributions)
K = 2;

% Number of times we will run the GMM (to avoid local minima)
n_runs = 100;

% Number of null folds
n_KLD = 50000;

% Tolerance threshold to stop the EM scheme
epsilon = 1e-5;

% Significance level (Bonferroni-corrected)
Alpha = 5/(n_regions+n_regions_SC);

warning('off');

% Actual computations for cortical regions
for r = 1:n_regions

    r
    
    % Our reference population is the set of HCs
    X_ref = squeeze(Morpho_HC(r,:,:));
    
    % Our other population is the set of baseline ET patients
    X = squeeze(Morpho_BASE(r,:,:));
     
    Result_actual_HB{r} = ET_DetermineSubcases(X,X_ref,K,n_runs,n_KLD,tol,Alpha);
end

% Same for non-cortical regions
for r = 1:n_regions_SC
    
    r
    
    % Our reference population is the set of HCs
    X_ref = squeeze(Morpho_HC_SC(r,:));
    
    % Our other population is the set of baseline ET patients
    X = squeeze(Morpho_BASE_SC(r,:));
    
    [Result_actual_SC_HB{r}] = ET_DetermineSubcases_1D(X,X_ref,K,n_runs,n_KLD,tol,Alpha);
end

% We sample the output data...
[GD_HB,LL_GMM_HB,LL_Gauss_HB,~,~,~,...
    ~,Acc_NB_HB,Acc_LDA_HB,dis_NB_HB,dis_LDA_HB,KLD_HB,KLDT_HB] = ...
    ET_GetDecision(Result_actual_HB,n_regions,1);

[GD_SC_HB,LL_GMM_SC_HB,LL_Gauss_SC_HB,~,~,~,~,Acc_NB_SC_HB,Acc_LDA_SC_HB,...
    dis_NB_SC_HB,~,KLD_SC_HB,KLDT_SC_HB] = ...
    ET_GetDecision(Result_actual_SC_HB,n_regions_SC,1);

% We extract the information regarding significant differences for the HC
% vs ET vs HC comparison
for r = 1:n_regions
    [pmat_HB(:,:,r),pvec_HB(:,r),p_KLD_HB(r)] = ET_GetRegionProfile(Result_actual_HB{r},Alpha);
end

for r = 1:n_regions_SC
    [pvec_HB_SC(:,r),p_KLD_HB_SC(r)] = ET_GetRegionProfile_1D(Result_actual_SC_HB{r});
end

for idx = 1:9
    idx_GD_HB(idx,:) = (pvec_HB(idx,:) < Alpha/100);
end

for idx = 1:2
    idx_GD_HB_SC(idx,:) = (pvec_HB_SC(idx,:) < Alpha/100);
end

idx_GD_HB = (sum(idx_GD_HB)>0);
idx_GD_HB_SC = (sum(idx_GD_HB_SC)>0);

% Regions for which the GMM wins over the multivariate Gaussian
idx_GMM_HB = find(LL_GMM_HB > LL_Gauss_HB & LL_GMM_HB < Inf);
idx_GMM_HB_SC = find(LL_GMM_SC_HB > LL_Gauss_SC_HB & LL_GMM_SC_HB < Inf);

if is_plot
    
    % Log-likelihood (cross-validated) comparison between Gaussian and GMM
    % cases
    figure;
    set(gca,'Box','off');
    bar(LL_Gauss_HB - LL_GMM_HB);
    
    figure;
    set(gca,'Box','off');
    bar(LL_Gauss_SC_HB - LL_GMM_SC_HB);
end

%%%%% The conclusion of this part is that the right putamen exhibit a 
%%%%% multimodal data structure. In addition, there
%%%%% are no group differences in terms of "mean" or covariance, but there 
%%%%% are in terms of variance



%% Part 5. Analysis of actual data (HC vs ET) in terms of log-likelihood, 
% to better understand the data organization

% Cortical regions...
for r = 1:n_regions

    % Data from the three groups
    tmp_data = squeeze(Morpho_HC(r,:,:));
    tmp_data2 = squeeze(Morpho_BASE(r,:,:));
    
    % Covariance and mean estimates in the HC population (maximum
    % likelihood estimates)
    COV(r,:,:) = cov(tmp_data);
    MU(r,:) = mean(tmp_data);
    
    % For control subjects, likelihood is computed for each region
    [L_HC(r,:),LL_HC(r,:)] = ET_EvaluateGaussian(tmp_data',MU(r,:)',squeeze(COV(r,:,:)));
    
    % Likelihood is also computed for the ET subjects, before and
    % after intervention, when assessing the probability that they could be
    % obtained from the HC distribution
    [L_BASE(r,:),LL_BASE(r,:)] = ET_EvaluateGaussian(tmp_data2',MU(r,:)',squeeze(COV(r,:,:)));
    
    % We perform the same by instead considering the pre-intervention ET
    % patients as reference
    COV2(r,:,:) = cov(tmp_data2);
    MU2(r,:) = mean(tmp_data2);
    
    % Evaluation of the baseline data wrt baseline distribution
    [L_tmp(r,:),LL_tmp(r,:)] = ET_EvaluateGaussian(tmp_data2',MU2(r,:)',squeeze(COV2(r,:,:)));
end

% Same on the non-cortical areas (only univariate)
for r = 1:n_regions_SC
    
    tmp_data = squeeze(Morpho_HC_SC(r,:));
    tmp_data2 = squeeze(Morpho_BASE_SC(r,:));

    COV_SC(r) = cov(tmp_data);
    MU_SC(r) = mean(tmp_data);
    
    [L_HC_SC(r,:),LL_HC_SC(r,:)] = ET_EvaluateGaussian(tmp_data,MU_SC(r),squeeze(COV_SC(r)));
    [L_BASE_SC(r,:),LL_BASE_SC(r,:)] = ET_EvaluateGaussian(tmp_data2,MU_SC(r),squeeze(COV_SC(r)));
   
    COV_SC2(r) = cov(tmp_data2);
    MU_SC2(r) = mean(tmp_data2);
    
    [L_tmp_SC(r,:),LL_tmp_SC(r,:)] = ET_EvaluateGaussian(tmp_data2,MU_SC2(r),COV_SC2(r));
end

% The more a subject is an "outlier" with respect to the dataset of ETpre
% values, the more that subject is impaired tremor-wise
X0 = [-LL_tmp',-LL_tmp_SC'];
Y0 = Symptoms_pre;
label = ones(n_BASE,1);
myPLS_inputs;
[IN_pre,pls_opts_pre] = myPLS_initialize(input,pls_opts,save_opts);
PLS_PRE_LL = myPLS_analysis(IN_pre,pls_opts_pre);

if is_plot
    
    figure;
    bar(PLS_PRE_LL.boot_results.Vb_mean(:,1));
    hold on
    errorbar(1:n_regions+n_regions_SC,PLS_PRE_LL.boot_results.Vb_mean(:,1),...
        PLS_PRE_LL.boot_results.Vb_mean(:,1)-PLS_PRE_LL.boot_results.Vb_lB(:,1),...
        PLS_PRE_LL.boot_results.Vb_uB(:,1)-PLS_PRE_LL.boot_results.Vb_mean(:,1),'color','k','LineStyle','None','CapSize',3);

    figure;
    bar(PLS_PRE_LL.boot_results.Ub_mean(:,1));
    hold on
    errorbar(1:size(Symptoms_pre,2),PLS_PRE_LL.boot_results.Ub_mean(:,1),...
        PLS_PRE_LL.boot_results.Ub_mean(:,1)-PLS_PRE_LL.boot_results.Ub_lB(:,1),...
        PLS_PRE_LL.boot_results.Ub_uB(:,1)-PLS_PRE_LL.boot_results.Ub_mean(:,1),'color','k','LineStyle','None');
end