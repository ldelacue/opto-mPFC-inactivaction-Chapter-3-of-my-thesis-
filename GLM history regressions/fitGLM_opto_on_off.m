function[Regressors,DependentVar, betas, pvalues, SEWeights, VIF] = fitGLM_opto_on_off(DataMatrix,LaserON_Ind,LaserOFF_Ind,OneHotEncoder,Dataframe,TrialsInPast)
%
% This function fits a Generalized Linear Model to the data formatted as
% DataMatrix (including withdrawals). The dependent variable (DependentVar)
% indicates R1 (0) vs R2 (1) response. The regressors used are 1 Current
% Stimulus Evidence 2a RewardHistoryBiases (t-1,t-2,t-3 & t-4) 2b
% UnrewardedHistoryBiases (t-1,t-2,t-3 & t-4) .

% GLM fitting: We asume that the data comes from a binomial distribution and
% the stretching function (also called link function or 1/link function) is logit
% function. This is also called a logistic regression.
% The regressors are coded from a R2 point of view. Meaning that a 1 means
% a positive tendency towards a R2 response and a -1  a negative one.

% Inputs:

% (1)   DataMatrix in allBDataKDB format
% (2)   Do we include the session by session regressor? 0 =no, 1= yes. If
%       you want to see the effect of the error trials you need 1=yes. This
%       is because in our task error influneces are so small that they you
%       need to clean the stronger session by session correlations that
%       there exist in our dataset so that you can actuall see those.
% (3)   Rats (Dataframe=0) or Pigeons (Dataframe = 2)?

% Book keeping:

% Luis de la Cuesta Ferrer   11/22. Script creation

if Dataframe ==2
    ValidTrials   = find(DataMatrix(:,3)~=0);    % For Pigeons
else
    ValidTrials   = not(isnan(DataMatrix(:,3))); % For Rats
end

[I,~] = find(DataMatrix(:,2)==DataMatrix(:,3));
DataMatrix(I,5)=0;

% Which indices you take: LASER ON TRIALS 

IndicesofInterest   = LaserON_Ind;

IndTrialsRewardedR1_ON_Log =  DataMatrix(IndicesofInterest,4)~=0 & DataMatrix(IndicesofInterest,3)==1  ;
IndTrialsRewardedR1_ON     =  IndicesofInterest(IndTrialsRewardedR1_ON_Log);
IndTrialsPunishedR1_ON_Log =  DataMatrix(IndicesofInterest,4)==0 & DataMatrix(IndicesofInterest,3)==1 ;
IndTrialsPunishedR1_ON     =  IndicesofInterest(IndTrialsPunishedR1_ON_Log);

IndTrialsRewardedR2_ON_Log =  DataMatrix(IndicesofInterest,4)~=0 & DataMatrix(IndicesofInterest,3)==2  ;
IndTrialsRewardedR2_ON     =  IndicesofInterest(IndTrialsRewardedR2_ON_Log);
IndTrialsPunishedR2_ON_Log =  DataMatrix(IndicesofInterest,4)==0 & DataMatrix(IndicesofInterest,3)==2 ;
IndTrialsPunishedR2_ON     =  IndicesofInterest(IndTrialsPunishedR2_ON_Log);

% Which indices you take: LASER OFF TRIALS 

IndicesofInterest   = LaserOFF_Ind;

IndTrialsRewardedR1_OFF_Log =  DataMatrix(IndicesofInterest,4)~=0 & DataMatrix(IndicesofInterest,3)==1  ;
IndTrialsRewardedR1_OFF     =  IndicesofInterest(IndTrialsRewardedR1_OFF_Log);
IndTrialsPunishedR1_OFF_Log =  DataMatrix(IndicesofInterest,4)==0 & DataMatrix(IndicesofInterest,3)==1 ;
IndTrialsPunishedR1_OFF     =  IndicesofInterest(IndTrialsPunishedR1_OFF_Log);

IndTrialsRewardedR2_OFF_Log =  DataMatrix(IndicesofInterest,4)~=0 & DataMatrix(IndicesofInterest,3)==2  ;
IndTrialsRewardedR2_OFF     =  IndicesofInterest(IndTrialsRewardedR2_OFF_Log);
IndTrialsPunishedR2_OFF_Log =  DataMatrix(IndicesofInterest,4)==0 & DataMatrix(IndicesofInterest,3)==2 ;
IndTrialsPunishedR2_OFF     =  IndicesofInterest(IndTrialsPunishedR2_OFF_Log);

for iTrialLag = 1:length(TrialsInPast)
    [LagTrialRewarded_ON(:,iTrialLag),LagTrialPunished_ON(:,iTrialLag)]  = BuildHistoryRegressors(iTrialLag,DataMatrix, IndTrialsRewardedR1_ON,IndTrialsRewardedR2_ON,IndTrialsPunishedR1_ON,IndTrialsPunishedR2_ON);
    [LagTrialRewarded_OFF(:,iTrialLag),LagTrialPunished_OFF(:,iTrialLag)]  = BuildHistoryRegressors(iTrialLag,DataMatrix, IndTrialsRewardedR1_OFF,IndTrialsRewardedR2_OFF,IndTrialsPunishedR1_OFF,IndTrialsPunishedR2_OFF);
end

% One hot encoding (session or condition wise).
% One categorical regressor per session that absorbs the session by session
% biases (it is 1 for that session trials and 0 for the rest)

if OneHotEncoder ==1 % Session
    BlockIdentifyer = zeros( length (DataMatrix), max(DataMatrix(:,8)) );
    for iSess = 1:max(DataMatrix(:,8))
        ThisSessionTrialIndices = find (DataMatrix(:,8)==iSess);
        BlockIdentifyer(ThisSessionTrialIndices,iSess) = 1;
    end
elseif OneHotEncoder ==2 % Condition
    BlockIdentifyer = zeros( length (DataMatrix), max(DataMatrix(:,6)) );
    for iCond = 1:max(DataMatrix(:,6))
        ThisConditionTrialIndices = find (DataMatrix(:,6)==iCond);
        BlockIdentifyer(ThisConditionTrialIndices,iCond) = 1;
    end
end
% Build regressors
% 1Current Stimulus Evidence (Dummy coded)

CurrentStimulus     = zeros(max(DataMatrix(:,1)),size(DataMatrix,1))';
for iStim = 1:max(DataMatrix(:,1))
    CurrentStimulus(DataMatrix(:,1)==iStim,iStim)=1;
end

% 2Y
DependentVar = (DataMatrix(:,3)-1);

% 3 Outcome history
if OneHotEncoder
    History_Regressors = [LagTrialRewarded_ON,LagTrialRewarded_OFF,LagTrialPunished_ON, LagTrialPunished_OFF, BlockIdentifyer(:,2:end)];
else
    History_Regressors = [LagTrialRewarded_ON,LagTrialRewarded_OFF,LagTrialPunished_ON, LagTrialPunished_OFF];

end

% 4a After effect module
% AfterEffect = zeros(1,size(DataMatrix,1));
% AfterEffect(2:end) = CurrentStimulus(1:end-1); % CurrentStimulus used to
% be the identity. Now it is a matrix dummy coding it.

% 5Regressors = [];

Regressors                  = [CurrentStimulus,History_Regressors];
Regressors                  = Regressors(ValidTrials,:);
DependentVar                = DependentVar(ValidTrials);


[b,dev,stats]               = glmfit(Regressors,DependentVar,'binomial','link','logit','Constant','off');
%[b,dev,stats]               = glmfit(Regressors,DependentVar,'binomial','link','logit');

betas     = stats.beta;
SEWeights = stats.se;
pvalues   = stats.p;

% Variance inflation coefficients (to check for colinearity of regressors)

%vif() computes variance inflation coefficients
%VIFs are also the diagonal elements of the inverse of the correlation matrix [1], a convenient result that eliminates the need to set up the various regressions
%[1] Belsley, D. A., E. Kuh, and R. E. Welsch. Regression Diagnostics. Hoboken, NJ: John Wiley & Sons, 1980.

R0 = corrcoef(Regressors); % correlation matrix
VIF  = diag(inv(R0))';
end


