%% Scatterplot_pR2_crit_RPs_vs_SPPs

% Compute change in biases (response probabilities and of criteria of the
% subjects) and plot them against each other to see how these two
% magnitudes relate. It additionally plots the change of PR2/PR1 and of
% criterion in Laser On vs OFF blocks to show the effect of the laser.

% Needs mSDTfit.m

%                                            Luis de la Cuesta 02/2025

clear all
close all

%% RPs
Epoch     = 1; % 1 During Stim 2 Outcome
Run_nr    = 1;
PR2_Laser_ON_All_RP  = [];
PR2_Laser_OFF_All_RP = [];
c_Laser_ON_All_RP    = [];
c_Laser_OFF_All_RP   = [];
Schedule_All_RP      = [];
Pheno_All_RP         = [];
SubjID_All_RP        = [];

Manipulation = 'RP';

if strcmp(Manipulation,'RP')
    load ("RP_Blockwise_opto_struct.mat")
    Condition = RP_Blockwise_opto_struct;
    BlockSize = 180;
    Punishment= 0;
    ProgrammedFirstTrialInSubBlock = [1,61,121, 181, 241, 301, 361, 421, 481, 541, 601, 661, 721, 781, 841]; % SubBlocks are datachunks of 60 trials. Normally, there are 3 subblocks per block. Each subblock will constitute a data point.

elseif strcmp(Manipulation,'SPP')
    load ("SPP_Blockwise_opto_struct.mat")
    Condition = SPP_Blockwise_opto_struct;
    BlockSize = 189;
    Punishment= 0;
    ProgrammedFirstTrialInSubBlock = [1,61,121, 181, 250, 310, 370, 439, 499, 559, 628, 688, 748, 817, 877]; % SubBlocks are datachunks of 60 trials. Normally, there are 3 subblocks per block. Each subblock will constitute a data point.

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
    [PR2_Laser_ON_Subj,PR2_Laser_OFF_Subj,PR2_Baseline_Laser_ON_Subj,PR2_Baseline_Laser_OFF_Subj,...
    c_Laser_ON_Subj,  c_Laser_OFF_Subj,  c_Baseline_Laser_ON_Subj,  c_Baseline_Laser_OFF_Subj]     = ComputeGlobalBiasMetrics (Datamatrix,Punishment);
    ID_Subj              = ones(size(PR2_Laser_ON_Subj,1),1)*iRow;
    Pheno_Subj           = ones(size(PR2_Laser_ON_Subj,1),1)*Pheno(iRow);


    PR2_Laser_ON_All_RP  = vertcat(PR2_Laser_ON_All_RP, PR2_Laser_ON_Subj-PR2_Baseline_Laser_ON_Subj);
    PR2_Laser_OFF_All_RP = vertcat(PR2_Laser_OFF_All_RP,PR2_Laser_OFF_Subj-PR2_Baseline_Laser_OFF_Subj);
    c_Laser_ON_All_RP    = vertcat(c_Laser_ON_All_RP, c_Laser_ON_Subj-c_Baseline_Laser_ON_Subj);
    c_Laser_OFF_All_RP   = vertcat(c_Laser_OFF_All_RP,c_Laser_OFF_Subj-c_Baseline_Laser_OFF_Subj);
    Schedule_All_RP      = vertcat(Schedule_All_RP,ones(size(PR2_Laser_ON_Subj,1),1)*Schedule(iRow));
    SubjID_All_RP =        vertcat(SubjID_All_RP,ID_Subj);
    Pheno_All_RP         = vertcat(Pheno_All_RP,Pheno_Subj);


    clearvars PR2_Laser_ON_Subj PR2_Laser_OFF_Subj PR2_Baseline_Laser_ON_Subj PR2_Baseline_Laser_OFF_Subj...
    c_Laser_ON_Subj  c_Laser_OFF_Subj  c_Baseline_Laser_ON_Subj  c_Baseline_Laser_OFF_Subj ID 
end
%% SPPs

Manipulation = 'SPP';
Epoch     = 1; % 1 During Stim 2 Outcome
Run_nr    = 1;
PR2_Laser_ON_All_SPP  = [];
PR2_Laser_OFF_All_SPP = [];
c_Laser_ON_All_SPP    = [];
c_Laser_OFF_All_SPP   = [];
Schedule_All_SPP      = [];
SubjID_All_SPP        = [];
Pheno_All_SPP         = [];


if strcmp(Manipulation,'RP')
    load ("RP_Blockwise_opto_struct.mat")
    Condition = RP_Blockwise_opto_struct;
    BlockSize = 180;
    Punishment= 0;
    ProgrammedFirstTrialInSubBlock = [1,61,121, 181, 241, 301, 361, 421, 481, 541, 601, 661, 721, 781, 841]; % SubBlocks are datachunks of 60 trials. Normally, there are 3 subblocks per block. Each subblock will constitute a data point.

elseif strcmp(Manipulation,'SPP')
    load ("SPP_Blockwise_opto_struct.mat")
    Condition = SPP_Blockwise_opto_struct;
    BlockSize = 189;
    Punishment= 0;
    ProgrammedFirstTrialInSubBlock = [1,61,121, 181, 250, 310, 370, 439, 499, 559, 628, 688, 748, 817, 877]; % SubBlocks are datachunks of 60 trials. Normally, there are 3 subblocks per block. Each subblock will constitute a data point.
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
    [PR2_Laser_ON_Subj,PR2_Laser_OFF_Subj,PR2_Baseline_Laser_ON_Subj,PR2_Baseline_Laser_OFF_Subj,...
     c_Laser_ON_Subj,  c_Laser_OFF_Subj,  c_Baseline_Laser_ON_Subj,  c_Baseline_Laser_OFF_Subj]     = ComputeGlobalBiasMetrics (Datamatrix,Punishment);
     ID_Subj              = ones(size(PR2_Laser_ON_Subj,1),1)*iRow;
     Pheno_Subj           = ones(size(PR2_Laser_ON_Subj,1),1)*Pheno(iRow);

    % Save individual data already Baseline normalized 
    PR2_Laser_ON_All_SPP  = vertcat(PR2_Laser_ON_All_SPP, PR2_Laser_ON_Subj-PR2_Baseline_Laser_ON_Subj);
    PR2_Laser_OFF_All_SPP = vertcat(PR2_Laser_OFF_All_SPP,PR2_Laser_OFF_Subj-PR2_Baseline_Laser_OFF_Subj);
    c_Laser_ON_All_SPP    = vertcat(c_Laser_ON_All_SPP, c_Laser_ON_Subj-c_Baseline_Laser_ON_Subj);
    c_Laser_OFF_All_SPP   = vertcat(c_Laser_OFF_All_SPP,c_Laser_OFF_Subj-c_Baseline_Laser_OFF_Subj);
    Schedule_All_SPP      = vertcat(Schedule_All_SPP,ones(size(PR2_Laser_ON_Subj,1),1)*Schedule(iRow));
    SubjID_All_SPP        = vertcat(SubjID_All_SPP,ID_Subj);
    Pheno_All_SPP         = vertcat(Pheno_All_SPP,Pheno_Subj);

    clearvars PR2_Laser_ON_Subj PR2_Laser_OFF_Subj PR2_Baseline_Laser_ON_Subj PR2_Baseline_Laser_OFF_Subj...
    c_Laser_ON_Subj  c_Laser_OFF_Subj  c_Baseline_Laser_ON_Subj  c_Baseline_Laser_OFF_Subj ID 
end

%% Plot
figure;

% Magnitude of bias elicited by RP and SPP manipulations
% swap inverted individuals (schedule ==2)

PR2_OFF_RP_inv                       = PR2_Laser_OFF_All_RP;
PR2_OFF_RP_inv(Schedule_All_RP==2)   = -1*PR2_Laser_OFF_All_RP(Schedule_All_RP==2);
PR2_ON_RP_inv                        = PR2_Laser_ON_All_RP;
PR2_ON_RP_inv(Schedule_All_RP==2)    = -1*PR2_Laser_ON_All_RP(Schedule_All_RP==2);

PR2_OFF_SPP_inv                      = PR2_Laser_OFF_All_SPP;
PR2_OFF_SPP_inv(Schedule_All_SPP==2) = -1*PR2_Laser_OFF_All_SPP(Schedule_All_SPP==2);
PR2_ON_SPP_inv                       = PR2_Laser_ON_All_SPP;
PR2_ON_SPP_inv(Schedule_All_SPP==2)  = -1*PR2_Laser_ON_All_SPP(Schedule_All_SPP==2);

c_OFF_RP_inv                         = c_Laser_OFF_All_RP;
c_OFF_RP_inv(Schedule_All_RP==1)     = -1*c_Laser_OFF_All_RP(Schedule_All_RP==1);
c_ON_RP_inv                          = c_Laser_ON_All_RP;
c_ON_RP_inv(Schedule_All_RP==1)      = -1*c_Laser_ON_All_RP(Schedule_All_RP==1);

c_OFF_SPP_inv                        = c_Laser_OFF_All_SPP;
c_OFF_SPP_inv(Schedule_All_SPP==1)   = -1*c_Laser_OFF_All_SPP(Schedule_All_SPP==1);
c_ON_SPP_inv                         = c_Laser_ON_All_SPP;
c_ON_SPP_inv(Schedule_All_SPP==1)    = -1*c_Laser_ON_All_SPP(Schedule_All_SPP==1);

r_RPs           = [ones(length(PR2_Laser_OFF_All_RP),1), PR2_Laser_OFF_All_RP]\c_Laser_OFF_All_RP;
ycalc_RP        = PR2_Laser_OFF_All_RP*r_RPs(2) +r_RPs(1);
scatter (PR2_Laser_OFF_All_RP,c_Laser_OFF_All_RP,'ok');hold on
model_RP        = fitlm(PR2_Laser_OFF_All_RP,c_Laser_OFF_All_RP);
p1              = plot(PR2_Laser_OFF_All_RP,ycalc_RP,'-k');hold on

r_SPPs          = [ones(length(PR2_Laser_OFF_All_SPP),1), PR2_Laser_OFF_All_SPP]\c_Laser_OFF_All_SPP;
ycalc_SPP       = PR2_Laser_OFF_All_SPP.*r_SPPs(2)+r_SPPs(1);
scatter (PR2_Laser_OFF_All_SPP,c_Laser_OFF_All_SPP,'ob');hold on
model_SPP       = fitlm(PR2_Laser_OFF_All_SPP,c_Laser_OFF_All_SPP);
p2              = plot(PR2_Laser_OFF_All_SPP,ycalc_SPP,'-b');hold on

axis square
axis([-0.5 0.5 -1.5 1.5])
ylabel('delta c')
xlabel('delta pR2')
title('Mapping c / pR2  - Laser OFF')
legend([p1,p2],{'RP manipulations','SPP manipulations'})


% Effect of Laser
%pR2 Halos
figure
title ('Halos')
ax1 = subplot(2,2,1);
scatter(PR2_OFF_RP_inv(Pheno_All_RP==1),PR2_ON_RP_inv(Pheno_All_RP==1),'x');hold on
line([-5,5],[-5,+5],'LineStyle', '--')
axis([-0.25 0.75 -0.25 0.75])
title (ax1,'delta PRbetter (RP) - HaloR')

axis square

ax2 = subplot(2,2,2);
scatter(c_OFF_RP_inv(Pheno_All_RP==1),c_ON_RP_inv(Pheno_All_RP==1),'x');hold on
line([-5,5],[-5,+5],'LineStyle', '--')
axis([-0.5 2.0 -0.5 2])
axis square
title (ax2, 'abs delta c (RP)')

%pR2 Halos

ax3 = subplot(2,2,3);
scatter(PR2_OFF_SPP_inv(Pheno_All_SPP==1),PR2_ON_SPP_inv(Pheno_All_SPP==1),'x');hold on
line([-5,5],[-5,+5],'LineStyle', '--')
axis([-0.25 0.75 -0.25 0.75])
xlabel('Laser OFF')
ylabel('Laser ON')
title (ax3, 'delta PRbetter (SPP)')


axis square
ax4 = subplot(2,2,4);
scatter(c_OFF_SPP_inv(Pheno_All_SPP==1),c_ON_SPP_inv(Pheno_All_SPP==1),'x');hold on
line([-5,5],[-5,+5],'LineStyle', '--')
axis([-0.5 2.0 -0.5 2])
axis square
title (ax4, 'abs delta c (SPP)')

% YFPs
figure
ax1 = subplot(2,2,1);
scatter(PR2_OFF_RP_inv(Pheno_All_RP==0),PR2_ON_RP_inv(Pheno_All_RP==0),'x');hold on
line([-5,5],[-5,+5],'LineStyle', '--')
axis([-0.25 0.75 -0.25 0.75])
title (ax1,'delta PRbetter (RP) - YFPs')
axis square

ax2 = subplot(2,2,2);
scatter(c_OFF_RP_inv(Pheno_All_RP==0),c_ON_RP_inv(Pheno_All_RP==0),'x');hold on
line([-5,5],[-5,+5],'LineStyle', '--')
axis([-0.5 2.0 -0.5 2])
axis square
title (ax2,'abs delta c (RP)')


ax3 = subplot(2,2,3);
scatter(PR2_OFF_SPP_inv(Pheno_All_SPP==0),PR2_ON_SPP_inv(Pheno_All_SPP==0),'x');hold on
line([-5,5],[-5,+5],'LineStyle', '--')
axis([-0.25 0.75 -0.25 0.75])
xlabel('Laser OFF')
ylabel('Laser ON')
title (ax3,'delta PRbetter (SPP)')


axis square
ax4 = subplot(2,2,4);
scatter(c_OFF_SPP_inv(Pheno_All_SPP==0),c_ON_SPP_inv(Pheno_All_SPP==0),'x');hold on
line([-5,5],[-5,+5],'LineStyle', '--')
axis([-0.5 2.0 -0.5 2])
axis square
title (ax4,'abs delta c (SPP)')


%% Stats

[h_RP_YFP_PR2_ON_OFF,p_RP_YFP_PR2_ON_OFF]     = ttest(PR2_OFF_RP_inv(Pheno_All_RP==0),PR2_ON_RP_inv(Pheno_All_RP==0));
[h_RP_YFP_c_ON_OFF,p_RP_YFP_c_ON_OFF]         = ttest(c_OFF_RP_inv(Pheno_All_RP==0),c_ON_RP_inv(Pheno_All_RP==0));
[h_RP_Halo_PR2_ON_OFF,p_RP_Halo_PR2_ON_OFF]   = ttest(PR2_OFF_RP_inv(Pheno_All_RP==1),PR2_ON_RP_inv(Pheno_All_RP==1));
[h_RP_Halo_c_ON_OFF,p_RP_Halo_c_ON_OFF]       = ttest(c_OFF_RP_inv(Pheno_All_RP==1),c_ON_RP_inv(Pheno_All_RP==1));
[h_SPP_Halo_PR2_ON_OFF,p_SPP_Halo_PR2_ON_OFF] = ttest(PR2_OFF_SPP_inv(Pheno_All_SPP==1),PR2_ON_SPP_inv(Pheno_All_SPP==1));
[h_SPP_Halo_c_ON_OFF,p_SPP_Halo_c_ON_OFF]     = ttest(c_OFF_SPP_inv(Pheno_All_SPP==1),c_ON_SPP_inv(Pheno_All_SPP==1));

% Levenes test (equal variance across criteria in SPPs vs RPs)
sample1 = ones(length(c_Laser_OFF_All_RP),1);
sample2 = 2*ones(length(c_Laser_OFF_All_SPP),1);
X_c = [c_Laser_OFF_All_RP,sample1;c_Laser_OFF_All_SPP,sample2];
alpha = 0.05;
Levenetest(X_c,alpha);
X_PR2 = [PR2_Laser_OFF_All_RP,sample1;PR2_Laser_OFF_All_SPP,sample2];
Levenetest(X_PR2,alpha);

%compare regression coefficients

x1 = PR2_Laser_OFF_All_RP;
y1 = c_Laser_OFF_All_RP;
x2 = PR2_Laser_OFF_All_SPP;
y2 = c_Laser_OFF_All_SPP;
[p_reg,stats] = mcompregslopes(x1,y1,x2,y2);



%% HELPER FUNCTIONS

%% function [] = ComputeGlobalBiasMetrics (Condition,RowsOfInterest)

function [PR2_Laser_ON_Subj,PR2_Laser_OFF_Subj,PR2_Baseline_Laser_ON_Subj,PR2_Baseline_Laser_OFF_Subj,...
          c_Laser_ON_Subj,  c_Laser_OFF_Subj,  c_Baseline_Laser_ON_Subj,  c_Baseline_Laser_OFF_Subj]     ...
          = ComputeGlobalBiasMetrics (Datamatrix,Punishment)

% Compute c and pR2 blockwise in Laser ON and OFF per subject. This function 
% is not really thought to be given different inputs just to do its thing and 
% not to hinder readability in the main program. Can be modified to return 
% different variables in the future (e.g: RT, MT, withdrawals...).                                

%                                               Luis 03/2025

Nr_subBlocks = 15;

    for iSess = 1:max(Datamatrix(:,18))
        SessionTrialInd              = find( Datamatrix(:,18) == iSess);
        SessionOptoTrialInd          = find( Datamatrix(SessionTrialInd,15)==1);
        if isempty(SessionOptoTrialInd) % There was not any stimulation in any block!
            Mssg =  "There was no laser stimulation in the session Nr"+ iSess;
            error(Mssg)
        end
        WichBlockIsStimulated(iSess) = Datamatrix(SessionOptoTrialInd(1),17);
    end
    clear iSess
    % Now lets compute some blockwise metrics of each session

    Alldprimes_Sess       = zeros(max(Datamatrix(:,18)),20);
    Allcriterions_Sess    = zeros(max(Datamatrix(:,18)),20);
    Allcriterions_1_Sess  = zeros(max(Datamatrix(:,18)),20);
    AllPR2_Sess           = zeros(max(Datamatrix(:,18)),20);
    AllResponseTimeR1     = zeros(max(Datamatrix(:,18)),20);
    AllwithdrawalRates    = zeros(max(Datamatrix(:,18)),20);


    for iSess = 1:max(Datamatrix(:,18)) % for all sessions

        ThisSessionTrials = find(Datamatrix(:,18)==iSess);
        ThisSessionData   = Datamatrix(ThisSessionTrials,:);
        LastTrial         = length(ThisSessionData);
        SessionData             = Datamatrix(Datamatrix(:,18)==iSess,:);
        [mu,c_from_model,r2]    = mSDTFit(SessionData,'single',60,'probit');
        Allcriterions_1_Sess(iSess,1:length(c_from_model)) = c_from_model;
        for jSubBlock = 1:max(Datamatrix(ThisSessionTrials,16)) % for all blocks in this session

            Data   = Datamatrix(Datamatrix(:,16) == jSubBlock & Datamatrix(:,18)==iSess,:);
            Rewards                =     length(Data(Data(:,9)==1));
            Withdrawals            =     length (find (isnan (Data (:,5))==1));
            CleanData              =     Data(not(isnan(Data(:,5))),:); %Clean NaN trials
            R1trials               =     CleanData(:,5)==1;
            R2trials               =     CleanData(:,5)==2;
            ResponseTimeR1         =     mean(CleanData(R1trials,8),1);
            ResponseTimeR2         =     mean(CleanData(R2trials,8),1);
            WithdrawalRate         =     Withdrawals/length(Data);
            RewardRate             =     Rewards/length(CleanData);
            Format                 =     1;
            ExtremeStimdprime      =     0;
            Categorywisedprime     =     0;
            [dprime,criterion,PR2] = compute_SDT_metrics(Data, Format,ExtremeStimdprime,Categorywisedprime );
            Alldprimes_Sess(iSess,jSubBlock)        = dprime;
            Allcriterions_Sess(iSess,jSubBlock)     = criterion;
            AllPR2_Sess (iSess,jSubBlock)           = PR2;
            AllResponseTimeR1(iSess, jSubBlock)       = ResponseTimeR1;
            AllResponseTimeR2(iSess, jSubBlock)       = ResponseTimeR2;
            AllWithdrawalRates_Sess (iSess,jSubBlock) = WithdrawalRate;
            AllRewardRates_Sess(iSess,jSubBlock)      = RewardRate;
        end


    end

    clear iSess jSubBlock dprime criterion PR2

    AllPR2_Sess             = AllPR2_Sess(:,1:Nr_subBlocks);
    Allcriterions_Sess      = Allcriterions_Sess(:,1:Nr_subBlocks);
    Allcriterions_1_Sess    = Allcriterions_1_Sess(:,1:Nr_subBlocks);
    Alldprimes_Sess         = Alldprimes_Sess(:,1:Nr_subBlocks);
    AllResponseTimeR1       = AllResponseTimeR1(:,1:Nr_subBlocks);
    AllResponseTimeR2       = AllResponseTimeR2(:,1:Nr_subBlocks);
    AllWithdrawalRates_Sess = AllWithdrawalRates_Sess(:,1:Nr_subBlocks);
    AllRewardRates_Sess     = AllRewardRates_Sess(:,1:Nr_subBlocks);

    if Punishment
        OptoInBlocks = [2,4];
    else
        OptoInBlocks = [3,5];
    end

    Meandprimes(1,:)         =  mean(Alldprimes_Sess((WichBlockIsStimulated==OptoInBlocks(1)),:),1);
    Meandprimes(2,:)         =  mean(Alldprimes_Sess((WichBlockIsStimulated==OptoInBlocks(2)),:),1);
    MeanPR2(1,:)             =  mean(AllPR2_Sess((WichBlockIsStimulated==OptoInBlocks(1)),:),1);
    MeanC(1,:)               =  mean(Allcriterions_Sess((WichBlockIsStimulated==OptoInBlocks(1)),:),1);
    MeanC_1(1,:)             =  mean(Allcriterions_1_Sess((WichBlockIsStimulated==OptoInBlocks(1)),:),1);
    PR2_Block3_ON            =  AllPR2_Sess(WichBlockIsStimulated==OptoInBlocks(1),:);
    PR2_Block5_ON            =  AllPR2_Sess(WichBlockIsStimulated==OptoInBlocks(2),:);
    c_Block3_ON              =  Allcriterions_1_Sess(WichBlockIsStimulated==OptoInBlocks(1),:);
    c_Block5_ON              =  Allcriterions_1_Sess(WichBlockIsStimulated==OptoInBlocks(2),:);
    MeanPR2(2,:)             =  mean(AllPR2_Sess((WichBlockIsStimulated==OptoInBlocks(2)),:),1);
    MeanC(2,:)               =  mean(Allcriterions_Sess((WichBlockIsStimulated==OptoInBlocks(2)),:),1);
    MeanC_1(2,:)             =  mean(Allcriterions_1_Sess((WichBlockIsStimulated==OptoInBlocks(2)),:),1);
    MeanResponseTimeR1(1,:)  =  mean(AllResponseTimeR1((WichBlockIsStimulated==OptoInBlocks(1)),:),1);
    MeanResponseTimeR1(2,:)  =  mean(AllResponseTimeR1((WichBlockIsStimulated==OptoInBlocks(2)),:),1);
    MeanResponseTimeR2(1,:)  =  mean(AllResponseTimeR2((WichBlockIsStimulated==OptoInBlocks(1)),:),1);
    MeanResponseTimeR2(2,:)  =  mean(AllResponseTimeR2((WichBlockIsStimulated==OptoInBlocks(2)),:),1);
    MeanWithdrawalRates(1,:) =  mean(AllWithdrawalRates_Sess((WichBlockIsStimulated==OptoInBlocks(1)),:),1);
    MeanWithdrawalRates(2,:) =  mean(AllWithdrawalRates_Sess((WichBlockIsStimulated==OptoInBlocks(2)),:),1);
    MeanRewardRates(1,:)     =  mean(AllRewardRates_Sess((WichBlockIsStimulated==OptoInBlocks(1)),:),1);
    MeanRewardRates(2,:)     =  mean(AllRewardRates_Sess((WichBlockIsStimulated==OptoInBlocks(2)),:),1);


    Block_1_datapoints = 1:1;
    Block_3_datapoints = 5:7;
    Block_5_datapoints = 11:13;

    PR2_Baseline_Laser_ON_Subj   = [mean(PR2_Block3_ON(:,Block_1_datapoints),2);mean(PR2_Block5_ON(:,Block_1_datapoints),2)];
    PR2_Baseline_Laser_OFF_Subj  = [mean(PR2_Block3_ON(:,Block_1_datapoints),2);mean(PR2_Block5_ON(:,Block_1_datapoints),2)];
    PR2_Laser_ON_Subj  = [mean(PR2_Block3_ON(:,Block_3_datapoints),2);mean(PR2_Block5_ON(:,Block_5_datapoints),2)];
    PR2_Laser_OFF_Subj = [mean(PR2_Block3_ON(:,Block_5_datapoints),2);mean(PR2_Block5_ON(:,Block_3_datapoints),2)];
    c_Baseline_Laser_ON_Subj  = [mean(c_Block3_ON(:,Block_1_datapoints),2);mean(c_Block5_ON(:,Block_1_datapoints),2)];
    c_Baseline_Laser_OFF_Subj = [mean(c_Block3_ON(:,Block_1_datapoints),2);mean(c_Block5_ON(:,Block_1_datapoints),2)];
    c_Laser_ON_Subj  = [mean(c_Block3_ON(:,Block_3_datapoints),2);mean(c_Block5_ON(:,Block_5_datapoints),2)];
    c_Laser_OFF_Subj = [mean(c_Block3_ON(:,Block_5_datapoints),2);mean(c_Block5_ON(:,Block_3_datapoints),2)];
 
end

%%
function [Levenetest] = Levenetest(X,alpha)
%Levene's Test for Homogeneity of Variances.
%[In the Levene's test the data are transforming to yij = abs[xij - mean(xj)]
%and uses the F distribution performing an one-way ANOVA using y as the 
%dependent variable (Brownlee, 1965; Miller, 1986)].
%
%   Syntax: function [Levenetest] = Levenetest(X,alfa) 
%      
%     Inputs:
%          X - data matrix (Size of matrix must be n-by-2; data=column 1, sample=column 2). 
%       alpha - significance level (default = 0.05).
%     Outputs:
%          - Sample variances vector.
%          - Whether or not the homoscedasticity was met.
%
%    Example: From the example 10.1 of Zar (1999, p.180), to test the Levene's
%             homoscedasticity of data with a significance level = 0.05.
%
%                                 Diet
%                   ---------------------------------
%                       1       2       3       4
%                   ---------------------------------
%                     60.8    68.7   102.6    87.9
%                     57.0    67.7   102.1    84.2
%                     65.0    74.0   100.2    83.1
%                     58.6    66.3    96.5    85.7
%                     61.7    69.8            90.3
%                   ---------------------------------
%                                       
%           Data matrix must be:
%            X=[60.8 1;57.0 1;65.0 1;58.6 1;61.7 1;68.7 2;67.7 2;74.0 2;66.3 2;69.8 2;
%            102.6 3;102.1 3;100.2 3;96.5 3;87.9 4;84.2 4;83.1 4;85.7 4;90.3 4];
%
%     Calling on Matlab the function: 
%             Levenetest(X)
%
%       Answer is:
%
% The number of samples are: 4
%
% ----------------------------
% Sample    Size      Variance
% ----------------------------
%   1        5         9.3920
%   2        5         8.5650
%   3        4         7.6567
%   4        5         8.3880
% ----------------------------
%   
% Levene's Test for Equality of Variances F=0.0335, df1= 3, df2=15
% Probability associated to the F statistic = 0.9914
% The associated probability for the F test is larger than 0.05
% So, the assumption of homoscedasticity was met.     
%

%  Created by A. Trujillo-Ortiz and R. Hernandez-Walls
%             Facultad de Ciencias Marinas
%             Universidad Autonoma de Baja California
%             Apdo. Postal 453
%             Ensenada, Baja California
%             Mexico.
%             atrujo@uabc.mx
%
%  April 11, 2003.
%
%  To cite this file, this would be an appropriate format:
%  Trujillo-Ortiz, A. and R. Hernandez-Walls. (2003). Levenetest: Levene's test for
%    homogeneity of variances. A MATLAB file. [WWW document]. URL http://www.mathworks.com/
%    matlabcentral/fileexchange/loadFile.do?objectId=3375&objectType=FILE
%
%  References:
% 
%  Brownlee, K. A. (1965) Statistical Theory and Methodology in Science
%           and Engineering. New York: John Wiley & Sons.
%  Miller, R. G. Jr. (1986) Beyond ANOVA, Basics of Applied Statistics.  
%           New York: John Wiley & Sons.
%  Zar, J. H. (1999), Biostatistical Analysis (2nd ed.).
%           NJ: Prentice-Hall, Englewood Cliffs. p. 180. 
%

if nargin < 2
   alpha = 0.05;
end 

Y=X;
k=max(Y(:,2));
fprintf('The number of samples are:%2i\n\n', k);

%Levene's Procedure.
n=[];s2=[];Z=[];
indice=Y(:,2);
for i=1:k
   Ye=find(indice==i);
   eval(['Y' num2str(i) '=Y(Ye,1);']);
   eval(['mY' num2str(i) '=mean(Y(Ye,1));']);
   eval(['n' num2str(i) '=length(Y' num2str(i) ') ;']);
   eval(['s2' num2str(i) '=(std(Y' num2str(i) ').^2) ;']);
   eval(['Z' num2str(i) '= abs((Y' num2str(i) ') - mY' num2str(i) ');']);
   eval(['xn= n' num2str(i) ';']);
   eval(['xs2= s2' num2str(i) ';']);
   eval(['x= Z' num2str(i) ';']);
   n=[n;xn];s2=[s2;xs2];Z=[Z;x];
end

Y=[Z Y(:,2)];

fprintf('-----------------------------\n');
disp(' Sample    Size      Variance')
fprintf('-----------------------------\n');
for i=1:k
   fprintf('   %d       %2i         %.4f\n',i,n(i),s2(i))
end
fprintf('-----------------------------\n');
disp(' ')

C=(sum(Y(:,1)))^2/length(Y(:,1)); %correction term.
SST=sum(Y(:,1).^2)-C; %total sum of squares.
dfT=length(Y(:,1))-1; %total degrees of freedom.

indice=Y(:,2);
for i=1:k
   Ye=find(indice==i);
   eval(['A' num2str(i) '=Y(Ye,1);']);
end

A=[];
for i=1:k
   eval(['x =((sum(A' num2str(i) ').^2)/length(A' num2str(i) '));']);
   A=[A,x];
end

SSA=sum(A)-C; %sample sum of squares.
dfA=k-1; %sample degrees of freedom.
SSE=SST-SSA; %error sum of squares.
dfE=dfT-dfA; %error degrees of freedom.
MSA=SSA/dfA; %sample mean squares.
MSE=SSE/dfE; %error mean squares.
F=MSA/MSE; %F-statistic.
v1=dfA;df1=v1;
v2=dfE;df2=v2;

P = 1 - fcdf(F,v1,v2);  %probability associated to the F-statistic.   

fprintf('Levene''s Test for Equality of Variances F=%3.4f, df1=%2i, df2=%2i\n', F,df1,df2);
fprintf('Probability associated to the F statistic = %3.4f\n', P);

if P >= alpha;
  fprintf('The associated probability for the F test is equal or larger than% 3.2f\n', alpha);
  fprintf('So, the assumption of homoscedasticity was met.\n');
else
  fprintf('The associated probability for the F test is smaller than% 3.2f\n', alpha);
  fprintf('So, the assumption of homoscedasticity was not met.\n');
end
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

%% Correct for extreme values (when there are 0 false alarms or 0 hits).
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

