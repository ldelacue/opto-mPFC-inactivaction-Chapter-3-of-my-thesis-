%% Linear Mixed Effect (LME) Model to check for the effect of mPFC inactivation on the absolute c
% Our dependent Variable,  c our fixed effect variable, Laser ON and Pheno; our random effect,  Rat ID.

%                                                       Luis de la Cuesta 02/2025

clear all
addpath '\\uni-mainz.de\dfs\Profiles\Settings\ldelacue\Desktop\code\Project 3 - Opto'
Condition = 1;
Epoch     = 1; % 1 During Stim 2 Outcome
Run_nr    = 1; % The RP experiment was run a second time (Run_nr = 2) only during stim, in 3 subjects (eNPHR) to show...
               % that the methodology was still effective 10 months later.
if Condition ==1
    load("RP_Blockwise_opto_struct.mat")
    Dataframe = RP_Blockwise_opto_struct;
    Punishement = 0;
    BlockSize = 180;
elseif Condition==2
    load("SPP_Blockwise_opto_struct.mat")
    Dataframe = SPP_Blockwise_opto_struct;
    Punishement = 0;
    BlockSize = 189;
elseif Condition==3
    load("PUN_Blockwise_opto_struct.mat")
    Dataframe = PUN_Blockwise_opto_struct;
    Punishement = 1;
end

All_AbsBias_vec     = [];
All_Bias_vec        = [];
All_Laser_vec       = [];
All_Phenotype_vec   = [];
All_Rat_vec         = [];
All_Sess_vec        = [];
All_OF_depth_vec    = [];

if     Epoch ==1
    RowsOfInterest = find( [Dataframe(:).LaserDur_ms] == 400 & [Dataframe(:).Run] == Run_nr); % DuringStim Inhibition 400ms

elseif Epoch==2
    RowsOfInterest = find( [Dataframe(:).LaserDur_ms] == 800 );                          % DuringOut  Inhibition 800ms
end

for iDataset=1:length(RowsOfInterest)
    iRow                   = RowsOfInterest(iDataset);
    allBDataKDB            = Dataframe(iRow).allBDataKDB;
    allBData               = Dataframe(iRow).allBData;
    SessRow                = Dataframe(iRow).allBData(:,18);
    Schedule               = Dataframe(iRow).BiasOrder;
    % Now lets compute some blockwise metrics of each session

    AllPR2_Sess           = zeros(1,max(SessRow));
    Alldprimes_Sess       = zeros(1,max(SessRow));
    Allcriterions_Sess    = zeros(1,max(SessRow));
    Allcriterions_1_Sess    = [];
    AllPR2_Sess           = zeros(1,max(SessRow));
    AllwithdrawalRates    = zeros(1,max(SessRow));


    for iSess = 1:max(SessRow) % for all sessions

        SessionData             = allBData(allBData(:,18)==iSess,:);
        [mu,c_from_model,r2]    = mSDTFit(SessionData,'single',60,'probit');
        Allcriterions_1_Sess(iSess,1:length(c_from_model)) = c_from_model;

        ThisSessionTrials = find(allBData(:,18)==iSess);
        Allcriterions_1_Sess = Allcriterions_1_Sess;
        for jSubBlock = 1:max(allBData(ThisSessionTrials,16)) % for all subblocks in this session

            Data   = allBData(allBData(:,16) == jSubBlock & allBData(:,18)==iSess,:);
            Withdrawals            =     length (find (isnan (Data (:,5))==1));
            WithdrawalRate         =     Withdrawals/length(Data);
            Format                 =     1;
            ExtremeStimdprime      =     1;
            Categorywisedprime     =     0;
            [dprime,criterion,PR2] = compute_SDT_metrics(Data, Format,ExtremeStimdprime,Categorywisedprime );
            Alldprimes_Sess(iSess,jSubBlock)        = dprime;
            Allcriterions_Sess(iSess,jSubBlock)     = criterion;
            AllPR2_Sess (iSess,jSubBlock)           = PR2;
            AllWithdrawalRates_Sess (iSess,jSubBlock)    = WithdrawalRate;
        end

    end

    clear iSess jSubBlock dprime criterion c_from_model PR2

    if Punishement
        SubBlocksofInterest = [6,12];
    else
        SubBlocksofInterest = [7,13]; % take directly PR2(or 1-PR2) at the blocks you are interested in comparing (namely 3 & 5).
    end


    if Punishement ~=1
        AbsBiasMatrix = Allcriterions_1_Sess(:,SubBlocksofInterest);
        if Schedule==1
            AbsBiasMatrix = AbsBiasMatrix*-1;
        elseif Schedule==2
            AbsBiasMatrix = AbsBiasMatrix;
        end
    elseif Punishement ==1
        if Schedule==1
            AbsBiasMatrix(:,1)   = Allcriterions_1_Sess(:,SubBlocksofInterest(1));
            AbsBiasMatrix(:,2)   = Allcriterions_1_Sess(:,SubBlocksofInterest(2))*-1;
        elseif Schedule==2
            AbsBiasMatrix(:,1)   = Allcriterions_1_Sess(:,SubBlocksofInterest(1))*-1;
            AbsBiasMatrix(:,2)   = Allcriterions_1_Sess(:,SubBlocksofInterest(2));
        end
    end

    % Laser
    LaserMatrix = zeros(max(SessRow),length(SubBlocksofInterest));

    for iSess = 1:max(allBData(:,18))
        SessionTrialInd              = find( allBData(:,18) == iSess);
        SessionOptoTrialInd          = find( allBData(SessionTrialInd,15)==1);
        if isempty(SessionOptoTrialInd) % There was not any stimulation in any block!
            Mssg =  "There was no laser stimulation in the session Nr"+ iSess;
            error(Mssg)
        end
        LaserMatrix_prime(iSess,allBData(SessionOptoTrialInd(1),17)) = 1; % was 1 instead of end
    end

    if Punishement
        LaserMatrix = LaserMatrix_prime(:,[2,4]);
    else
        LaserMatrix = LaserMatrix_prime(:,[3,5]);
    end
    clear iSess

    % Rat & Pheno & Optic_F_depth

    OF_depth_matrix = ones(max(SessRow),length(SubBlocksofInterest));
    OF_depth = (Dataframe(iRow).DV_OpticFiber_Planned);
    OF_depth_matrix = OF_depth_matrix*OF_depth;

    RatIDMatrix = ones(max(SessRow),length(SubBlocksofInterest))*iDataset;
    if strcmp(Dataframe(iRow).Pheno,"Halo")
        PhenotypeMatrix  = ones(max(SessRow),length(SubBlocksofInterest));
    else
        PhenotypeMatrix  = zeros(max(SessRow),length(SubBlocksofInterest));
    end

    % Session &
    SessionIDMatrix = ones(max(SessRow),length(SubBlocksofInterest));
    for iSess = 1:max(allBData(:,18))
        SessionIDMatrix(iSess,:) = SessionIDMatrix(iSess,:)*iSess;
    end

    %% Visualize the data you aim to model
    % figure
    % Laser3SessInd =find(LaserMatrix(:,1)==1);
    % Laser5SessInd =find(LaserMatrix(:,2)==1);
    % labels = ({'Early + Laser','Late','Early','Late + Laser'});
    % boxplot([AbsBiasMatrix(Laser3SessInd,:),AbsBiasMatrix(Laser5SessInd,:)],labels)
    % title("Phenotype " + Dataframe(iRow).Pheno )
    % PlotOptoBlockData(allBData,Punishement, BlockSize, BiasOrder)
    %%
    All_AbsBias_vec     = [All_AbsBias_vec,AbsBiasMatrix(:)'];
    All_Laser_vec       = [All_Laser_vec,LaserMatrix(:)'];
    All_Phenotype_vec   = [All_Phenotype_vec,PhenotypeMatrix(:)'];
    All_Rat_vec         = [All_Rat_vec,RatIDMatrix(:)'];
    All_Sess_vec        = [All_Sess_vec, SessionIDMatrix(:)'];
    All_OF_depth_vec    = [All_OF_depth_vec, OF_depth_matrix(:)'];
    clearvars AbsBiasMatrix LaserMatrix PhenotypeMatrix RatIDMatrix SessionIDMatrix OF_depth_matrix LaserMatrix_prime

end

% Generate Table
VariableNames= {'AbsoluteBias','laser','phenotype', 'session','rat','OF_depth'};
tbl = table( All_AbsBias_vec', All_Laser_vec', All_Phenotype_vec', categorical(All_Sess_vec'), categorical(All_Rat_vec'), All_OF_depth_vec', 'VariableNames',VariableNames);
%% Test hypothesis
%% 1) First Run: Does the block with laser have an effect in pR2 in Controls

if Run_nr ==1
    % % Fit Model
    lme = fitlme(tbl,'AbsoluteBias ~ 1 + laser + phenotype + laser*phenotype + (1|rat)')

    % PLot results
    PR2    = table2array(tbl(:,1));
    laser  = table2array(tbl(:,2));
    pheno  = table2array(tbl(:,3));
    session = double(table2array(tbl(:,4)));
    ratID   = double(table2array(tbl(:,5)));
    nrats   = numel(unique(ratID));
    figure
    colors = lines(nrats);
    subplot(231),hold on
    boxplot(PR2,laser)
    for i = 1:nrats,for j = 1:2
            plot(j+0.3+rand*0.1,PR2(ratID==i & laser==j-1),'.','Color',colors(i,:))
    end,end
axis([.5,2.5,-0.5,1.5]),ylabel('criterion')
set(gca,'XTick',1:2,'XTickLabel',{'Laser off','Laser on'})
subplot(232),hold on
boxplot(PR2,pheno)
for j = 1:2
    plot(j+0.3+rand*0.1,PR2(pheno==j-1),'.','Color',colors(i,:))
end
set(gca,'XTick',1:2,'XTickLabel',{'eYFP','eNpHR'})
subplot(233),hold on
boxplot(PR2,[pheno,laser])
for i = 1:nrats
    try
        plot(0.95+rand*0.1,PR2(ratID==i & pheno==0 & laser==0),'k.')
    catch
    end
    try
        plot(1.95+rand*0.1,PR2(ratID==i & pheno==0 & laser==1),'k.')
    catch
    end
    try
        plot(2.95+rand*0.1,PR2(ratID==i & pheno==1 & laser==0),'k.')
    catch
    end
    try
        plot(3.95+rand*0.1,PR2(ratID==i & pheno==1 & laser==1),'k.')
    catch
    end
    hold on

end
set(gca,'XTick',1:4,'XTickLabel',{'eYFP, laser off','eYFP, laser on','eNpHR, laser off','eNpHR, laser on'})
subplot(234),title('Mean and SEM'),hold on
c = 1;
for i = 1:2
    for j = 1:2
        errorbar(c,mean(PR2(pheno==i-1 & laser==j-1)),sem(PR2(pheno==i-1 & laser==j-1)),'.','CapSize',0)
        c = c+1;
    end
end
set(gca,'XTick',1:4,'XTickLabel',{'eYFP, laser off','eYFP, laser on','eNpHR, laser off','eNpHR, laser on'})
axis([.5,4.5,-0.25,0.75])
ylabel('criterion')
subplot(235),title('Median and median AD'),hold on
c = 1;
for i = 1:2
    for j = 1:2
        errorbar(c,median(PR2(pheno==i-1 & laser==j-1)),mad(PR2(pheno==i-1 & laser==j-1),1),'.','CapSize',0)
        c = c+1;
    end
end
set(gca,'XTick',1:4,'XTickLabel',{'eYFP, laser off','eYFP, laser on','eNpHR, laser off','eNpHR, laser on'})
axis([.5,4.5,-0.25,0.75])

%% 2) Second Run, no Controls
elseif Run_nr ==2
    lme = fitlme(tbl,'AbsoluteBias ~ 1 + laser + (1|rat)')

    % Plot results
    PR2    = table2array(tbl(:,1));
    laser  = table2array(tbl(:,2));
    session = double(table2array(tbl(:,4)));
    ratID   = double(table2array(tbl(:,5)));
    nrats   = numel(unique(ratID));
    close all
    figure
    colors = lines(nrats);
    subplot(131);
    hold on
    boxplot(PR2,laser)
    for i = 1:nrats,for j = 1:2
            plot(j-0.05+rand*0.1,PR2(ratID==i & laser==j-1),'.','Color',colors(i,:))
    end,end
axis([.5,4.5,-0.5,1.5]),ylabel('criterion')
set(gca,'XTick',1:2,'XTickLabel',{'Laser off','Laser on'})

subplot(132),title('Mean and SEM'),hold on
c = 1;
for j = 1:2
    errorbar(c,mean(PR2(laser==j-1)),sem(PR2(laser==j-1)),'.','CapSize',0)
    c = c+1;
end
set(gca,'XTick',1:2,'XTickLabel',{'laser off','laser on'})
axis([.5,4.5,-0.5,1.5])
ylabel('criterion')
subplot(133),title('Median and median AD'),hold on
c = 1;
for j = 1:2
    errorbar(c,median(PR2(laser==j-1)),mad(PR2(laser==j-1),1),'.','CapSize',0)
    c = c+1;
end

set(gca,'XTick',1:2,'XTickLabel',{'laser off','laser on'})
axis([.5,4.5,-0.5,1.5])


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


