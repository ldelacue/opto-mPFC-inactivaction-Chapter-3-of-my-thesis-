%% BlockTrials_Psych_Curve_pbetter_allSubjects

% PLots the mean psychommetric curve on Baseline blocks vs exp. manipulation
% blocks when Laser is ON vs OFF for eNpHR/eYFP- expressing rats, for the
% condition of interest (Manipulation ='RP','SPP' or 'PUN');

%                                                       Luis de la Cuesta 02/2025

clear all
Manipulation = 'SPP';
Epoch     = 1; % 1 During Stim 2 Outcome
Run_nr    = 1;

if strcmp(Manipulation,'RP')
    load ("RP_Blockwise_opto_struct.mat")
    Condition = RP_Blockwise_opto_struct;

elseif strcmp(Manipulation,'SPP')
    load ("SPP_Blockwise_opto_struct.mat")
    Condition = SPP_Blockwise_opto_struct;

elseif strcmp(Manipulation,'PUN')
    load ("Pun_Blockwise_opto_struct.mat")
    Condition = PUN_Blockwise_opto_struct;
end

if     Epoch ==1
    RowsOfInterest = find( [Condition(:).LaserDur_ms] == 400 & [Condition(:).Run] == Run_nr); % DuringStim Inhibition 400ms
elseif Epoch==2
    RowsOfInterest = find( [Condition(:).LaserDur_ms] == 800 );                          % DuringOut  Inhibition 800ms
end

for iRow=1:length(RowsOfInterest)

    Schedule(iRow) =Condition(RowsOfInterest(iRow)).BiasOrder;          % 1 or 2
    Datamatrix = Condition(RowsOfInterest(iRow)).allBData;
    Pheno(iRow)      = strcmp(Condition(RowsOfInterest(iRow)).Pheno,'Halo');
    % Which block has stimulation
    for iSess = 1:max(Datamatrix(:,18))
        SessionTrialInd              = find( Datamatrix(:,18) == iSess);
        SessionOptoTrialInd          = find( Datamatrix(SessionTrialInd,15)==1);
        if isempty(SessionOptoTrialInd) % There was not any stimulation in any block!
            Mssg =  "There was no laser stimulation in the session Nr"+ iSess;
            error(Mssg)
        end
        WichBlockIsStimulated(iSess) = Datamatrix(SessionOptoTrialInd(1),17);
        LaserMatrix(iSess,Datamatrix(SessionOptoTrialInd(1),17)) = 1;
    end
    clear iSess
    %% Now lets compute some stimulus-wise Pr2 for each block of interest ( baseline 60 trials, manipulation 120 trials)

contingencies = {'Baseline','Bias','Bias Opto'} ;
UpallBData                 = Datamatrix(not(isnan(Datamatrix(:,5))),:);
for iCont = 1:length(contingencies)
switch contingencies{iCont}
case 'Baseline'
    SubBlockOfInterest         = 1; 
    TrialsOfInterest = find(UpallBData(:,16)==SubBlockOfInterest);
    [Leftw_baseline,Outcome_baseline] = PrepPsychommetricTable (  UpallBData, TrialsOfInterest);
    pR2_baseline(iRow,:) = Leftw_baseline(4,:);
    plot(Leftw_baseline(4,:))
    axis([0.5 6.5 0 1])
    axis square
    case 'Bias'
    
    SubBlockOfInterest        = [6,7,12,13]; 
    AllControlTrials = UpallBData(UpallBData(:,15)==0,:);
    TrialsOfInterest = find(AllControlTrials(:,16)==SubBlockOfInterest(1)|AllControlTrials(:,16)==SubBlockOfInterest(2)|...
                            AllControlTrials(:,16)==SubBlockOfInterest(3)|AllControlTrials(:,16)==SubBlockOfInterest(4));
    [Leftw_bias_control,Outcome_bias_control] = PrepPsychommetricTable (  AllControlTrials, TrialsOfInterest);
    pR2_bias_control(iRow,:) = Leftw_bias_control(4,:);
    hold on; plot(Leftw_bias_control(4,:))
case 'Bias Opto'
    
    SubBlockOfInterest        = [6,7,12,13]; 
    AllOptoTrials = UpallBData(UpallBData(:,15)==1,:);
    TrialsOfInterest = find(AllOptoTrials(:,16)==SubBlockOfInterest(1)|AllOptoTrials(:,16)==SubBlockOfInterest(2)|...
                            AllOptoTrials(:,16)==SubBlockOfInterest(3)|AllOptoTrials(:,16)==SubBlockOfInterest(4));
    [Leftw_bias_opto,Outcome_bias_opto] = PrepPsychommetricTable (  AllOptoTrials, TrialsOfInterest);
    pR2_bias_opto(iRow,:) = Leftw_bias_opto(4,:);
    hold on;plot(Leftw_bias_opto(4,:))
end
end


    %% Premature withdrawals
    LabelConditions = 0;

    %
  
  clearvars  Alldprimes_Sess Allcriterions_Sess AllPR2_Sess AllResponseTimeR1 AllResponseTimeR2 AllWithdrawalRates_Sess Meandprimes ...
      MeanPR2 WichBlockIsStimulated
end

pR_better_baseline(Schedule==1,:) = pR2_baseline(Schedule==1,:);
pR_better_baseline(Schedule==2,:) = flip( 1-pR2_baseline(Schedule==2,:), 2);
pR_better_opto(Schedule==1,:)     = pR2_bias_opto(Schedule==1,:);
pR_better_opto(Schedule==2,:)     = flip( 1-pR2_bias_opto(Schedule==2,:), 2);
pR_better_control(Schedule==1,:)  = pR2_bias_control(Schedule==1,:);
pR_better_control(Schedule==2,:)  = flip( 1-pR2_bias_control(Schedule==2,:), 2);


sem_baseline_Halos = std(pR_better_baseline(Pheno==1,:),1)./sqrt(length(pR_better_baseline(Pheno==1)));
sem_baseline_YFPs  = std(pR_better_control(Pheno==0,:),1)./sqrt(length(pR_better_control(Pheno==0)));
sem_control_Halos  = std(pR_better_control(Pheno==1,:),1)./sqrt(length(pR_better_control(Pheno==1)));
sem_control_YFPs   = std(pR_better_control(Pheno==0,:),1)./sqrt(length(pR_better_control(Pheno==0)));
sem_opto_Halos     = std(pR_better_opto(Pheno==1,:),1)./sqrt(length(pR_better_opto(Pheno==1)));
sem_opto_YFPs      = std(pR_better_opto(Pheno==0,:),1)./sqrt(length(pR_better_opto(Pheno==0)));

%% Plotting aggregated data
orange= [0 0.4470 0.7410];
blue=   [0.8500 0.3250 0.0980];
figure
if strcmp(Manipulation,'SPP')
StimulusValues =  Condition(4).StimSequence;
else
StimulusValues =  Condition(1).StimSequence;
end
subplot(1,2,1)

errorbar(StimulusValues,mean(pR_better_baseline(Pheno==1,:)),sem_baseline_Halos,'-','LineWidth',2); hold on
errorbar(StimulusValues,mean(pR_better_control(Pheno==1,:)),sem_control_Halos,'-k','LineWidth',2); hold on
errorbar(StimulusValues,mean(pR_better_opto(Pheno==1,:)),sem_opto_Halos,'-','color',[0.9290 0.6940 0.1250],'LineWidth',2)
axis([0 100 0 1])
axis square
ylabel('pR better')
xlabel('Morphs')
legend('Baseline','Bias Laser OFF',' Bias Laser ON')
title('P. curves (Halo) n=',num2str(length(find(Pheno==1))) )

subplot(1,2,2)
errorbar(StimulusValues,mean(pR_better_baseline(Pheno==0,:)),sem_baseline_YFPs,'-','LineWidth',2); hold on
errorbar(StimulusValues,mean(pR_better_control(Pheno==0,:)),sem_control_YFPs,'-k','LineWidth',2); hold on
errorbar(StimulusValues,mean(pR_better_opto(Pheno==0,:)),sem_opto_YFPs,'-','color',[0.9290 0.6940 0.1250],'LineWidth',2)
axis([0 100 0 1])
axis square
xlabel('Morphs')
title('P. curves (YFP) n=',num2str(length(find(Pheno==0))) )


%% HELPER FUNCTIONS

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
