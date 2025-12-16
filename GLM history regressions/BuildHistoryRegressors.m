 function[LastTrialRewarded,LastTrialPunished] = BuildHistoryRegressors(TrialsInPast,allBDataKDB, IndTrialsRewardedR1,IndTrialsRewardedR2,IndTrialsPunishedR1,IndTrialsPunishedR2)

 % This function builds the regressors of the influence of a given
 % outcome positive and negative for the response
 % specified at trial = TrialInPast to investigate its effect on 
 % the current decisionin the current decision. It needs as input
 % IndTrialsRewardedR1/R2 and IndTrialsPunishedR1/R2 .
 
 %                                                  Luis 21/10/2022
 


    LastTrialRewarded = zeros(1,size(allBDataKDB,1)+TrialsInPast);
    LastTrialPunished = zeros(1,size(allBDataKDB,1)+TrialsInPast);

    LastTrialRewarded(1,IndTrialsRewardedR1+TrialsInPast) =-1; %Response 1 Rewarded trials            %  A note on the regressors signs: If you give the same sign to R1 correct and R1 incorrect
    LastTrialRewarded(1,IndTrialsRewardedR2+TrialsInPast) =1; %Response 2 Rewarded trials                you would expect negative regressors for the unrewarded outcomes. This visually makes 
                                                                                                      %  more sense but clashes with our negative learning rates interpretation of the error-based
    LastTrialPunished(1,IndTrialsPunishedR1+TrialsInPast) =-1; %Response 1 Punished trials             %  models.
    LastTrialPunished(1,IndTrialsPunishedR2+TrialsInPast) =1; %Response 2 Punished trials

    LastTrialRewarded= LastTrialRewarded(:,1:length(LastTrialRewarded)-TrialsInPast)'; %reformatting of the matrix
    LastTrialPunished= LastTrialPunished(:,1:length(LastTrialPunished)-TrialsInPast)'; %reformatting of the matrix
    



 
