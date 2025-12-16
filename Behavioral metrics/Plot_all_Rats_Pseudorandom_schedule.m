%% Plot_all_Rats_Pseudorandom_schedule
% Calculate general behaviour metrics for trials which feature Laser ON vs
% OFF and perform mean comparison. In hind sight, I believe a LME should be
% run instead of pulling all trials together per subject but did not get to
% do it. 

% Needs msep (Maik's function)
                                           % Luis de la Cuesta 02/2025

clear all
load('Pseudorandom_opto_struct.mat')
RowsOfInterest_Halo =  [];
RowsOfInterest_YFP  =  [];
Exclude_first_subjects =1; 
Implement_drift_correction = 1; % Implement criterion drift correction (1) or not (0) as in Lak et al., 2020

% We excluded the first subjects because they underwent testing on the pseudorandom condition 
% after having gone through the blocked RP manipulations sessions
% and we found that animals seemed to have associated optic stimulation
% with the response associated to the contingency with which it was
% previously associated in the previous experiment (RP condition). We decided to not interpret
% the data of these initial animals for these reason. In the rest of
% animals, they were allways tested before in the pseudorandom condition
% prior to the blocked RP manipulations. 

if     Exclude_first_subjects ==1
PseudorandomOptoStruct = PseudorandomOptoStruct_beta(:,4:end);
else
    PseudorandomOptoStruct = PseudorandomOptoStruct_beta;
end

for iRow=1:length(PseudorandomOptoStruct)
    if strcmp(PseudorandomOptoStruct(iRow).Pheno,'Halo')
      RowsOfInterest_Halo =  [RowsOfInterest_Halo, iRow];
    else strcmp(PseudorandomOptoStruct(iRow).Pheno,'YFP')
      RowsOfInterest_YFP =  [RowsOfInterest_YFP, iRow];
    end
end


RowsOfInterest = RowsOfInterest_YFP;
for iRat=1:length(RowsOfInterest)
%% Compute Probs of potential interest

DataMatrix = PseudorandomOptoStruct(RowsOfInterest(iRat)).allBData;
StimulusWitdrawalsInd = DataMatrix(:,13)==2|DataMatrix(:,13)==4; % Withdrawals in which the stimuli already had an identity
CompletedTrialsInd    = not(isnan(DataMatrix(:,5)));

for iStim = 1:max(DataMatrix(:,3))
    WithSOptoTrials(iStim,iRat) = length(DataMatrix(StimulusWitdrawalsInd & DataMatrix(:,3)==iStim & DataMatrix(:,15)==1));
    WithSControlTrials(iStim,iRat) = length(DataMatrix(StimulusWitdrawalsInd & DataMatrix(:,3)==iStim & DataMatrix(:,15)==0));
    CompletedSOptoTrials(iStim,iRat) = length(DataMatrix(CompletedTrialsInd & DataMatrix(:,3)==iStim & DataMatrix(:,15)==1));
    CompletedSControlTrials(iStim,iRat) = length(DataMatrix(CompletedTrialsInd & DataMatrix(:,3)==iStim & DataMatrix(:,15)==0));
end

WithROptoTrials(1,iRat) = length(DataMatrix(StimulusWitdrawalsInd & DataMatrix(:,15)==1))/length(DataMatrix(DataMatrix(:,15)==1));
WithRControlTrials(1,iRat) = length(DataMatrix(StimulusWitdrawalsInd & DataMatrix(:,15)==0))/length(DataMatrix(DataMatrix(:,15)==0));
WithROptoTrials(2,iRat) = length(DataMatrix(DataMatrix(:,13)==2 & DataMatrix(:,15)==1)) /length(DataMatrix(StimulusWitdrawalsInd & DataMatrix(:,15)==1)); % Proportion of S1 withdrawals for withdrawals with stimulus identity
WithRControlTrials(2,iRat) = length(DataMatrix(DataMatrix(:,13)==2 & DataMatrix(:,15)==0)) /length(DataMatrix(StimulusWitdrawalsInd & DataMatrix(:,15)==0)); % Proportion of S2 withdrawals for withdrawals with stimulus identity
WithROptoTrials(3,iRat) = length(DataMatrix(DataMatrix(:,13)==4 & DataMatrix(:,15)==1)) /length(DataMatrix(StimulusWitdrawalsInd & DataMatrix(:,15)==1)); % Proportion of S1 withdrawals for withdrawals with stimulus identity
WithRControlTrials(3,iRat) = length(DataMatrix(DataMatrix(:,13)==4 & DataMatrix(:,15)==0)) /length(DataMatrix(StimulusWitdrawalsInd & DataMatrix(:,15)==0)); % Proportion of S2 withdrawals for withdrawals with stimulus identity

DataMatrix = DataMatrix (not(isnan(DataMatrix(:,5))),:); % Get rid of withdrawals
RTs = DataMatrix(:,7)-DataMatrix(:,4)-DataMatrix(:,8);
RT(1,iRat) = mean(RTs(DataMatrix(:,15)==0));
RT(2,iRat) = mean(RTs(DataMatrix(:,15)==1));
MT(1,iRat) = mean(DataMatrix(DataMatrix(:,15)==0,8));
MT(2,iRat) = mean(DataMatrix(DataMatrix(:,15)==1,8));

StimulusValues(iRat,:) = PseudorandomOptoStruct(RowsOfInterest(iRat)).StimSequence;
%SquaredStimulusValues  = StimulusValues(iRat,:).^2;
%ScaledSquaredSimulusValues(iRat,:) = SquaredStimulusValues./max(SquaredStimulusValues)*100;
%StimulusValues(iRat,:) = ScaledSquaredSimulusValues(iRat,:);

% Correct

OptoTrials                = find(DataMatrix(:,15)==1);
ControlTrials             = find(DataMatrix(:,15)==0);

RewardedTrials            = find(DataMatrix(:,9)==1 );
UnrewardedTrials          = find(DataMatrix(:,9)~=1 );
RewardsR1                 = find(DataMatrix(:,9)==1 & DataMatrix(:,5)==1);
RewardsR2                 = find(DataMatrix(:,9)==1 & DataMatrix(:,5)==2);
Before_Rew_Right             = RewardsR1(2:end)-1;
Before_Rew_Left             = RewardsR2(2:end)-1;

% dprime and crit
ExtremeStimdprime=1;
Categorywisedprime=0;
Format =1;
Data = DataMatrix(ControlTrials,:);
[dprime(1,iRat),criterion(1,iRat),PR2(1,iRat),Cat1HR(1,iRat),Cat2FAR(1,iRat)] = compute_SDT_metrics(Data,Format,ExtremeStimdprime,Categorywisedprime);
Data = DataMatrix(OptoTrials,:);
[dprime(2,iRat),criterion(2,iRat),PR2(2,iRat),Cat1HR(2,iRat),Cat2FAR(2,iRat)] = compute_SDT_metrics(Data,Format,ExtremeStimdprime,Categorywisedprime);




pROptoTrials(1,iRat) = length(DataMatrix(DataMatrix(:,5)==1 & DataMatrix(:,15)==1))/length(DataMatrix(DataMatrix(:,15)==1));
pRControlTrials(1,iRat) = length(DataMatrix(DataMatrix(:,5)==1 & DataMatrix(:,15)==0))/length(DataMatrix(DataMatrix(:,15)==0));
pROptoTrials(2,iRat) = length(DataMatrix(DataMatrix(:,5)==2 & DataMatrix(:,15)==1))/length(DataMatrix(DataMatrix(:,15)==1));
pRControlTrials(2,iRat) = length(DataMatrix(DataMatrix(:,5)==2 & DataMatrix(:,15)==0))/length(DataMatrix(DataMatrix(:,15)==0));

HROptoTrials(1,iRat) = length(DataMatrix(DataMatrix(:,6)==1 & DataMatrix(:,15)==1))/length(DataMatrix(DataMatrix(:,15)==1));
HRControlTrials(1,iRat) = length(DataMatrix(DataMatrix(:,6)==1 & DataMatrix(:,15)==0))/length(DataMatrix(DataMatrix(:,15)==0));


RewardedOptoTrials       = find(DataMatrix(:,9)==1 & DataMatrix(:,15)==1);
UnrewardedOptoTrials     = find(DataMatrix(:,9)~=1 & DataMatrix(:,15)==1);

RewardedOptoTrialsRight       = find(DataMatrix(:,5)==1 & DataMatrix(:,9)==1 & DataMatrix(:,15)==1);
RewardedOptoTrialsLeft       = find(DataMatrix(:,5)==2 & DataMatrix(:,9)==1 & DataMatrix(:,15)==1);

RewardedControlTrialsRight       = find(DataMatrix(:,5)==1 & DataMatrix(:,9)==1 & DataMatrix(:,15)==0);
RewardedControlTrialsLeft       = find(DataMatrix(:,5)==2 & DataMatrix(:,9)==1 & DataMatrix(:,15)==0);

RewardedControlTrials    = find(DataMatrix(:,9)==1 & DataMatrix(:,15)==0);
UnrewardedControlTrials  = find(DataMatrix(:,9)~=1 & DataMatrix(:,15)==0);

PercRewardedOpto      = length(RewardedOptoTrials)/length(OptoTrials);
PercUnrewardedControl = length(UnrewardedControlTrials)/length(ControlTrials);

% Compute leftward probabilities for different type of trials

AfterRewardedTrials       = 1+RewardedTrials (1:end-1);
AfterUnrewardedTrials     = 1+UnrewardedTrials (1:end-1);

AfterOptoTrials           = OptoTrials(1:end-1)+1;
BeforeOptoTrials          = OptoTrials(2:end)-1;

FirstHalfOptoTrials       = OptoTrials(1:ceil( length(OptoTrials)/2 ));
SecondHalfOptoTrials      = OptoTrials(ceil( length(OptoTrials)/2 ): length(OptoTrials));

AfterOptoRewardedTrials    = 1+RewardedOptoTrials(1:end-1);
AfterOptoUnrewardedTrials  = 1+UnrewardedOptoTrials(1:end-1);

AfterOptoRewardedRightTrials   = 1+RewardedOptoTrialsRight(1:end-1);
AfterOptoRewardedLeftTrials    = 1+RewardedOptoTrialsLeft(1:end-1);

AfterControlRewardedRightTrials   = 1+RewardedControlTrialsRight(1:end-1);
AfterControlRewardedLeftTrials    = 1+RewardedControlTrialsLeft(1:end-1);

OptoAfterRewardedTrial     = 1+find(DataMatrix(BeforeOptoTrials,9)==1);
OptoAfterUnrewardedTrial   = 1+find(DataMatrix(BeforeOptoTrials,9)~=1);

AfterControlTrials        = ControlTrials(1:end-1)+1;
BeforeControlTrials       = ControlTrials(2:end)-1;

FirstHalfControlTrials    = ControlTrials(1:ceil( length(ControlTrials)/2 ));
SecondHalfControlTrials   = ControlTrials(ceil( length(ControlTrials)/2 ): length(ControlTrials));

AfterControlRewardedTrials    = 1+RewardedControlTrials(1:end-1);
AfterControlUnrewardedTrials  = 1+UnrewardedControlTrials(1:end-1);

    %ALL TRIALS

    TrialsOfInterest          = 1:size(DataMatrix,1);
    [AllTrialsLeftw,~]        = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
     AllTrials_pLeftw(iRat,:) = AllTrialsLeftw(4,:);
    %BEFORE R1 REW TRIALS (to remove effect of drift)
    TrialsOfInterest = Before_Rew_Right;
    [BeforeRightRewLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
     BeforeRightRew_pLeftw(iRat,:) = BeforeRightRewLeftw(4,:);


    %BEFORE R1 REW TRIALS (to remove effect of drift)
    TrialsOfInterest = Before_Rew_Left;
    [BeforeLeftRewLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
     BeforeLeftRew_pLeftw(iRat,:) = BeforeLeftRewLeftw(4,:);

    %AFTER REWARDED TRIALS
    
    TrialsOfInterest = AfterRewardedTrials;
    [AfterRewardedLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);

    %AFTER UNREWARDED TRIALS
    
    TrialsOfInterest = AfterUnrewardedTrials;
    [AfterUnrewardedLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);

    %OPTO

    TrialsOfInterest = OptoTrials;
    [OptoLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
    OptoPLeftw(iRat,:) = OptoLeftw(4,:);

    %TRIAL AFTER OPTO

    TrialsOfInterest = AfterOptoTrials;
    [AfterOptoLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
    
  
    % TRIAL BEFORE OPTO

    TrialsOfInterest = BeforeOptoTrials;
    [BeforeOptoLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);

    % 1st HALF OPTO TRIALS
    
    TrialsOfInterest = FirstHalfOptoTrials;
    [FirstHalfOptoLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);

    % 2nd HALF OPTO TRIALS

    TrialsOfInterest = SecondHalfOptoTrials;
    [SecondHalfOptoLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
    
    % AFTER OPTO REWARDED TRIAL

    TrialsOfInterest = AfterOptoRewardedTrials;
    [AfterOptoRewardedLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);

    % AFTER OPTO REWARDED RIGHT TRIAL

    TrialsOfInterest = AfterOptoRewardedRightTrials;
    [AfterOptoRewardedRight_Leftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
     AfterOptoRewardedRight_pLeftw(iRat,:) = AfterOptoRewardedRight_Leftw(4,:);

    % AFTER OPTO REWARDED LEFT TRIAL

    TrialsOfInterest = AfterOptoRewardedLeftTrials;
    [AfterOptoRewardedLeft_Leftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
     AfterOptoRewardedLeft_pLeftw(iRat,:) = AfterOptoRewardedLeft_Leftw(4,:);

    % AFTER CONTROL REWARDED RIGHT TRIAL

    TrialsOfInterest = AfterControlRewardedRightTrials;
    [AfterControlRewardedRight_Leftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
     AfterControlRewardedRight_pLeftw(iRat,:) = AfterControlRewardedRight_Leftw(4,:);

    % AFTER CONTROL REWARDED LEFT TRIAL

    TrialsOfInterest = AfterControlRewardedLeftTrials;
    [AfterControlRewardedLeft_Leftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
     AfterControlRewardedLeft_pLeftw(iRat,:) = AfterControlRewardedLeft_Leftw(4,:);

    % AFTER OPTO UNREWARDED TRIAL

    TrialsOfInterest = AfterOptoUnrewardedTrials;
    [AfterOptoUnrewardedLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);

    % CONTROL (NOT OPTO)
    
    TrialsOfInterest = ControlTrials;
    [ControlLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
    ControlPLeftw(iRat,:) = ControlLeftw(4,:);
  
    %TRIAL AFTER CONTROL (NOT OPTO)

    TrialsOfInterest = AfterControlTrials;
    [AfterControlLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
    
    %TRIAL BEFORE CONTROL (NOT OPTO)

    TrialsOfInterest = BeforeControlTrials;
    [BeforeControlLeftw, ~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);

    % 1st HALF CONTROL TRIALS

    TrialsOfInterest = FirstHalfControlTrials;
    [FirstHalfControlLeftw, ~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);

    % 2nd HALF CONTROL TRIALS

    TrialsOfInterest = SecondHalfControlTrials;
    [SecondHalfControlLeftw, ~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);

    % AFTER CONTROL REWARDED TRIAL

    TrialsOfInterest = AfterControlRewardedTrials;
    [AfterControlRewardedLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);

    % AFTER CONTROL UNREWARDED TRIAL

    TrialsOfInterest = AfterControlUnrewardedTrials;
    [AfterControlUnrewardedLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
    
    % OPTO AFTER REWARDED TRIAL
     TrialsOfInterest =  OptoAfterRewardedTrial;
    [OptoAfterRewardedLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);
   
    
    % OPTO AFTER UNREWARDED TRIAL
    
     TrialsOfInterest =  OptoAfterUnrewardedTrial;
    [OptoAfterUnrewardedLeftw,~]  = PrepPsychommetricTable ( DataMatrix, TrialsOfInterest);

Data = DataMatrix(AfterControlRewardedLeftTrials,:);
[dprime_after_rew_Left(1,iRat),criterion_after_rew_Left(1,iRat),~,~,~] = compute_SDT_metrics(Data,Format,ExtremeStimdprime,Categorywisedprime);
Data = DataMatrix(AfterOptoRewardedLeftTrials,:);
[dprime_after_rew_Left(2,iRat),criterion_after_rew_Left(2,iRat),~,~,~] = compute_SDT_metrics(Data,Format,ExtremeStimdprime,Categorywisedprime);

Data = DataMatrix(AfterControlRewardedRightTrials,:);
[dprime_after_rew_Right(1,iRat),criterion_after_rew_Right(1,iRat),~,~,~] = compute_SDT_metrics(Data,Format,ExtremeStimdprime,Categorywisedprime);
Data = DataMatrix(AfterOptoRewardedRightTrials,:);
[dprime_after_rew_Right(2,iRat),criterion_after_rew_Right(2,iRat),~,~,~] = compute_SDT_metrics(Data,Format,ExtremeStimdprime,Categorywisedprime);
end

%% PLot
figure; 
subplot(4,2,1)

ControlPLeftw(end+1,:) = mean(ControlPLeftw,1);
ControlPLeftw(end+1,:) = std(ControlPLeftw(1:end-1,:),1)/sqrt(size(ControlPLeftw,1)-1); % 4 is the number of datapoints 

OptoPLeftw(end+1,:)    = mean(OptoPLeftw,1);
OptoPLeftw(end+1,:)    = std(OptoPLeftw(1:end-1,:),1)/sqrt(size(OptoPLeftw,1)-1); % 4 is the number of datapoints 


p1=plot   (StimulusValues(1,:),ControlPLeftw(end-1,:),'k-','LineWidth',2.5);hold all ;
errorbar(StimulusValues(1,:),ControlPLeftw(end-1,:),ControlPLeftw(end,:),'k.','MarkerSize',20)
%Opto
p2=plot   (StimulusValues(1,:),OptoPLeftw(end-1,:),'y-','LineWidth',2.5);hold all ;
errorbar(StimulusValues(1,:),OptoPLeftw(end-1,:),OptoPLeftw(end,:),'y.','MarkerSize',20)

xlabel('Morph (defined as fraction High Freq chord)')
ylabel('Fraction leftward')
axis([0 100 0 1])
axis square

title('Effect on current trial')
legend([p1,p2],{'Laser OFF','Laser ON'})


% Plot PR2 contingent on rewarded choice on previous trial for Opto and Control conditions
if Implement_drift_correction
  DriftBeforeLeftRew  = abs(AllTrials_pLeftw-BeforeLeftRew_pLeftw);
  DriftBeforeRightRew = abs(AllTrials_pLeftw-BeforeRightRew_pLeftw);
  AfterControlRewardedLeft_pLeftw  = AfterControlRewardedLeft_pLeftw  - DriftBeforeLeftRew; 
  AfterOptoRewardedLeft_pLeftw     = AfterOptoRewardedLeft_pLeftw     - DriftBeforeLeftRew; 
  AfterControlRewardedRight_pLeftw = AfterControlRewardedRight_pLeftw + DriftBeforeRightRew; 
  AfterOptoRewardedRight_pLeftw    = AfterOptoRewardedRight_pLeftw    + DriftBeforeRightRew; 
end
subplot(4,2,2)

AfterControlRewardedLeft_pLeftw(end+1,:)  = mean(AfterControlRewardedLeft_pLeftw,1);
AfterControlRewardedLeft_pLeftw(end+1,:)  = std(AfterControlRewardedLeft_pLeftw(1:end-1,:),1)/sqrt(size(AfterControlRewardedLeft_pLeftw,1)-1); % 4 is the number of datapoints 
AfterControlRewardedRight_pLeftw(end+1,:) = mean(AfterControlRewardedRight_pLeftw,1);
AfterControlRewardedRight_pLeftw(end+1,:) = std(AfterControlRewardedRight_pLeftw(1:end-1,:),1)/sqrt(size(AfterControlRewardedRight_pLeftw,1)-1); % 4 is the number of datapoints 

AfterOptoRewardedLeft_pLeftw(end+1,:)  = mean(AfterOptoRewardedLeft_pLeftw,1);
AfterOptoRewardedLeft_pLeftw(end+1,:)  = std(AfterOptoRewardedLeft_pLeftw(1:end-1,:),1)/sqrt(size(AfterOptoRewardedLeft_pLeftw,1)-1); % 4 is the number of datapoints 
AfterOptoRewardedRight_pLeftw(end+1,:) = mean(AfterOptoRewardedRight_pLeftw,1);
AfterOptoRewardedRight_pLeftw(end+1,:) = std(AfterOptoRewardedRight_pLeftw(1:end-1,:),1)/sqrt(size(AfterOptoRewardedRight_pLeftw,1)-1); % 4 is the number of datapoints 

plot   (StimulusValues(1,:),AfterControlRewardedLeft_pLeftw(end-1,:),'k--','LineWidth',2.5);hold all ;
errorbar(StimulusValues(1,:),AfterControlRewardedLeft_pLeftw(end-1,:),AfterControlRewardedLeft_pLeftw(end,:),'k.','MarkerSize',20)
plot   (StimulusValues(1,:),AfterControlRewardedRight_pLeftw(end-1,:),'k:','LineWidth',2.5);hold all ;
errorbar(StimulusValues(1,:),AfterControlRewardedRight_pLeftw(end-1,:),AfterControlRewardedRight_pLeftw(end,:),'k.','MarkerSize',20)

plot   (StimulusValues(1,:),AfterOptoRewardedLeft_pLeftw(end-1,:),'y--','LineWidth',2.5);hold all ;
errorbar(StimulusValues(1,:),AfterOptoRewardedLeft_pLeftw(end-1,:),AfterOptoRewardedLeft_pLeftw(end,:),'y.','MarkerSize',20)
plot   (StimulusValues(1,:),AfterOptoRewardedRight_pLeftw(end-1,:),'y:','LineWidth',2.5);hold all ;
errorbar(StimulusValues(1,:),AfterOptoRewardedRight_pLeftw(end-1,:),AfterOptoRewardedRight_pLeftw(end,:),'y.','MarkerSize',20)

xlabel('Morph (defined as fraction High Freq chord)')
ylabel('Fraction leftward')
axis([0 100 0 1])
axis square

title('Effect on subsequent trial')
legend('After leftw. reward','','After rightw. reward')

%PLot criterion and dprime

subplot(4,2,3)
plot(criterion);hold on; plot(mean(criterion,2),'k-*')
axis square
title('c')
xticklabels({'Laser OFF','','Laser ON'})
axis([0.5,2.5,-0.5,0.5])

subplot(4,2,4)
plot(dprime); hold on; plot(mean(dprime,2),'k-*')
axis square
title('dprime')
xticklabels({'Laser OFF','','Laser ON'})
axis([0.5,2.5,0,3])

% PLot Completed & Withdrawn Trials, Stimulus-wise (no difference)

subplot(4,2,5) % Controls
plot(WithSControlTrials,'k--'); hold all; plot(CompletedSControlTrials,'b--');
p1=plot(mean(WithSControlTrials,2),'k-','LineWidth',3); p2=plot(mean(CompletedSControlTrials,2),'b-','LineWidth',3);
xlabel('Stim Nr')
ylabel('Trials')
legend([p1,p2],{'Nr. withdrawn trials with identitiy','Nr. completed trials'})
title('Control Trials')
axis([0.5 6.5 0 600])

axis square
subplot(4,2,6) % Opto
plot(WithSOptoTrials,'k--'); hold all; plot(CompletedSOptoTrials,'b--');
p1=plot(mean(WithSOptoTrials,2),'k-','LineWidth',3); p2=plot(mean(CompletedSOptoTrials,2),'b-','LineWidth',3);
xlabel('Stim Nr')
ylabel('Trials')
legend([p1,p2],{'Nr. withdrawn trials with identitiy','Nr. completed trials'})
title('Opto Trials')
axis([0.5 6.5 0 600])
axis square

% PLot MT & RT

subplot(4,2,7)
plot(MT);hold on;plot(mean(MT,2),'k-*')
axis square
title('Movement times')
xticklabels({'Laser OFF','','Laser ON'})
axis([0.5,2.5,0.3,0.9])

subplot(4,2,8)
plot(RT);hold on;plot(mean(RT,2),'k-*')
axis square
title('Reaction times')
xticklabels({'Laser OFF','','Laser ON'})
axis([0.5,2.5,0.16,0.22])



%% Stats


[h_RT,p_RT] = ttest(RT(1,:),RT(2,:));
[h_MT,p_MT] = ttest(MT(1,:),MT(2,:));

[h_c,p_c] = ttest(criterion(1,:),criterion(2,:));
[h_d,p_d] = ttest(dprime(1,:),dprime(2,:));

[h_c_after_rew_L,p_c_after_rew_L] = ttest(criterion_after_rew_Left(1,:),criterion_after_rew_Left(2,:));   %R2
[h_c_after_rew_R,p_c_after_rew_R] = ttest(criterion_after_rew_Right(1,:),criterion_after_rew_Right(2,:)); %R1

[h_d_after_rew_L,p_d_after_rew_L] = ttest(dprime_after_rew_Left(1,:),dprime_after_rew_Left(2,:));
[h_d_after_rew_R,p_d_after_rew_R] = ttest(dprime_after_rew_Right(1,:),dprime_after_rew_Right(2,:));

%% HELPER FUNCTIONS
%% function [Leftw,Outcome] = PrepPsychommetricTable ( UpallBData, TrialsOfInterest)

function [Leftw,Outcome] = PrepPsychommetricTable ( UpallBData, TrialsOfInterest)

for iStim=1:max(UpallBData(:,3)) % Loop over stim

    Leftw(1,iStim)        = length(find(UpallBData(TrialsOfInterest,3)==iStim & UpallBData(TrialsOfInterest,5)==2));     %find nr of RII (leftwards) responses for each stim
    Leftw(2,iStim)        = length(find(UpallBData(TrialsOfInterest,3)==iStim & UpallBData(TrialsOfInterest,5)==1));     %find nr of RI (rightwards) responses for each stim
    
    Outcome(1,iStim)      = length(find(UpallBData(TrialsOfInterest,3)==iStim & UpallBData(TrialsOfInterest,6)==1));     %find nr of hits for each stim
    Outcome(2,iStim)      = length(find(UpallBData(TrialsOfInterest,3)==iStim & UpallBData(TrialsOfInterest,6)==0));     %find nr of FA for each stim

        % Correction for cases with exclusive response. Create half a response
    % for the unchosen response. This avoids infinits values in perceptual
    % space calculations.
    
    if     Leftw(1,iStim) == 0
        Leftw(1,iStim) = 0.5;

    elseif Leftw(2,iStim) == 0
        Leftw(2,iStim) = 0.5;

    end

end

Leftw(3,:)=sum(Leftw,1);
Leftw(4,:)=Leftw(1,:)./Leftw(3,:);
Leftw(5,:)=msep(Leftw(4,:),Leftw(3,:));
end

%% function[dprime,criterion,PR2,Cat1HR,Cat2FAR] = compute_SDT_metrics(Data,Format,ExtremeStimdprime,Categorywisedprime)

function[dprime,criterion,PR2,Cat1HR,Cat2FAR] = compute_SDT_metrics(Data,Format,ExtremeStimdprime,Categorywisedprime)

% Inputs: Data can be provided in allBDataKDB (Format =0) or allBData
% (Format =1) forms. Unless otherwise specified (Categorywisedprime ==1) it defines the category of the
% stimuli assuming that there are 6 stimuli, 1-3 are Cat 1 and 4-6 are Cat 2. Returns d prime and c.
% It returns only one value so you need to feed the single datapoint for
% which you want the SDT metrics
                                                                            % Luis de la Cuesta  09/22
% ExtremeStimdprime: If 1, dprime is calculated only using stim 1 and stim
% 6 (the easiest stimuli). Otherwise, it encompasses the performance in all
% 6 stimuli.
% Categorywisedprime =1. Use column 2 to define what categories trials are
% instead of the Stimuli. This is preferable by default
                                                                            % Luis de la Cuesta  10/23

if Categorywisedprime ==1
    ExtremeStimdprime = 0;
end

if Format ==0

    Data = Data(not(isnan(Data(:,3))),:); %Clean NaN trials in Rats
    Data = Data(~Data(:,3)==0,:);         %Clean NaN trials in Rats. These two do not interphere
    if Categorywisedprime ==1

        Cat1TrialsInd        = find(Data(:,2)==1);
        Cat2TrialsInd        = find(Data(:,2)==2);
    elseif Categorywisedprime ==0
        if ExtremeStimdprime ==0
            Cat1TrialsInd        = find(Data(:,1)==1|Data(:,1)==2|Data(:,1)==3  );
            Cat2TrialsInd        = find(Data(:,1)==4|Data(:,1)==5|Data(:,1)==6  );
        elseif ExtremeStimdprime ==1
            ExtremeStim = [min(Data(:,1)) max(Data(:,1))];
            Cat1TrialsInd        = find(Data(:,1)==ExtremeStim(1) );
            Cat2TrialsInd        = find(Data(:,1)==ExtremeStim(2) );
        end
    end
    Cat1Hits             = length( find (Data(Cat1TrialsInd,2)==Data(Cat1TrialsInd,3) ));
    Cat1FA               = length( find (Data(Cat1TrialsInd,2)~=Data(Cat1TrialsInd,3) ));

    Cat2Hits             = length( find (Data(Cat2TrialsInd,2)==Data(Cat2TrialsInd,3) ));
    Cat2FA               = length( find (Data(Cat2TrialsInd,2)~=Data(Cat2TrialsInd,3) ));
    
    NrR2                  = length( find (Data(:,3)==2 ) );



elseif Format ==1

    %Data = Data(min(Data(:,1)):max(Data(:,1)),:);      % Remove trials after last responded trial
    Data = Data(not(isnan(Data(:,5))),:); % Clean NaN trials

    Cat1TrialsInd        = find(Data(:,3)==1|Data(:,3)==2|Data(:,3)==3  );
    Cat2TrialsInd        = find(Data(:,3)==4|Data(:,3)==5|Data(:,3)==6  );

    Cat1Hits             = length( find (Data(Cat1TrialsInd,5)==1 ));
    Cat1FA               = length( find (Data(Cat1TrialsInd,5)==2 ));

    Cat2Hits             = length( find (Data(Cat2TrialsInd,5)==2 ));
    Cat2FA               = length( find (Data(Cat2TrialsInd,5)==1 ));

    NrR2                  = length( find (Data(:,5)==2 ) );

    if ExtremeStimdprime
        Cat1TrialsInd        = find(Data(:,3)==1 );
        Cat2TrialsInd        = find(Data(:,3)==6 );

        Cat1Hits             = length( find (Data(Cat1TrialsInd,5)==1 ));
        Cat1FA               = length( find (Data(Cat1TrialsInd,5)==2 ));

        Cat2Hits             = length( find (Data(Cat2TrialsInd,5)==2 ));
        Cat2FA               = length( find (Data(Cat2TrialsInd,5)==1 ));
    end


end


Cat1HR               = Cat1Hits/length(Cat1TrialsInd);
Cat1FAR              = Cat1FA/length(Cat1TrialsInd);

Cat2HR               = Cat2Hits/length(Cat2TrialsInd);
Cat2FAR              = Cat2FA/length(Cat2TrialsInd);

PR2                  = NrR2 / length (Data);

% Correct for extreme values (when there are 0 false alarms or 0 hits).
% Instead use 1/2N. Half a trial to avoid inifity values. Read more at:

% Brown and White 2005.
% Behaviour Research methods


if Cat1Hits == 0
    Cat1HR  = 0.5 *(1/length(Cat1TrialsInd));
    Cat1FAR = 1- 0.5 *(1/length(Cat1TrialsInd));

elseif Cat1FA == 0
    Cat1HR = 1- 0.5 *(1/length(Cat1TrialsInd));
    Cat1FAR = 0.5 *(1/length(Cat1TrialsInd));
end

if Cat2Hits == 0
    Cat2HR = 0.5 *(1/length(Cat2TrialsInd));
    Cat2FAR = 1- 0.5 *(1/length(Cat2TrialsInd));

elseif Cat2FA == 0
    Cat2HR = 1- 0.5 *(1/length(Cat2TrialsInd));
    Cat2FAR = 0.5 *(1/length(Cat2TrialsInd));
end

dprime               = norminv(Cat1HR)-norminv(Cat2FAR);
criterion            = 0.5*(norminv(Cat1HR)+ norminv(Cat2FAR));

end

