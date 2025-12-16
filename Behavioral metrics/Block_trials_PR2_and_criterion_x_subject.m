%% Block_trials PR2 and criterion x subject
% PLot the data PR2 and criterion for a given subject (animalID). PLots
% both RP and PUN data. PLots both mPFC inactivation during stimulus (Epoch=1) and
% outcome (Epoch=1). The cirterion being plotted follows linear (plotted first) and probit 
% (plotted second) fits.

%                                                       Luis de la Cuesta 02/2025

clear all
close all
%% load dataset (individual subject, exp. manipulation and opto ephoc)
Manipulation = 'RP';
Epoch     = 1; % 1 During Stim 2 Outcome
animalID  =  'LF75';

if strcmp(Manipulation,'RP')
    load ("RP_Blockwise_opto_struct.mat")
    Condition = RP_Blockwise_opto_struct;
    BlockSize = 180;
    Punishement= 0;

elseif strcmp(Manipulation,'SPP')
    load ("SPP_Blockwise_opto_struct.mat")
    Condition = SPP_Blockwise_opto_struct;
    BlockSize = 189;
    Punishement= 0;

elseif strcmp(Manipulation,'PUN')
    load ("Pun_Blockwise_opto_struct.mat")
    Condition = PUN_Blockwise_opto_struct;
    BlockSize = 240;
    Punishement= 1;
end

if Punishement==1
    BlockChanges         = [120, 360, 480, 720, 1720 ];                  % (global variables: these variables do not change)
    BlockLengths         = [120, 240, 120, 240, 1000  ];
    ProgrammedFirstTrialInBlock    = [1,  121,  361, 481, 721 ];                  % Blocks are session divisions in which  an experimental variable is manipulated (in this case we manipulate RP1 and RP2 blockwise)
    ProgrammedFirstTrialInSubBlock = [1,61,121, 181, 241, 301, 361, 421, 481, 541, 601, 661,721,781]; % SubBlocks are datachunks of 60 trials. Normally, there are 3 subblocks per block. Each subblock will constitute a data point.

elseif BlockSize == 189 % For SPP

    BlockChanges         = [60, 60+BlockSize*1, 60+BlockSize*2, 60+BlockSize*3, 60+BlockSize*4, 1780 ];                  % (global variables: these variables do not change)
    BlockLengths         = [60, 189, 189, 189, 189, 960  ];
    ProgrammedFirstTrialInBlock    = [1,  BlockChanges(1)+1,  BlockChanges(2)+1, BlockChanges(3)+1, BlockChanges(4)+1, 781  ];                  % Blocks are session divisions in which  an experimental variable is manipulated (in this case we manipulate RP1 and RP2 blockwise)
    ProgrammedFirstTrialInSubBlock = [1,61,121, 181, 250, 310, 370, 439, 499, 559, 628, 688, 748, 817, 877]; % SubBlocks are datachunks of 60 trials. Normally, there are 3 subblocks per block. Each subblock will constitute a data point.

elseif BlockSize == 180 % For RP

    BlockChanges         = [60, 240, 420, 600, 780, 1780 ];                  % (global variables: these variables do not change)
    BlockLengths         = [60, 180, 180, 180, 180, 960  ];
    ProgrammedFirstTrialInBlock    = [1,  61,  241, 421, 601, 781  ];                  % Blocks are session divisions in which  an experimental variable is manipulated (in this case we manipulate RP1 and RP2 blockwise)
    ProgrammedFirstTrialInSubBlock = [1,61,121, 181, 241, 301, 361, 421, 481, 541, 601, 661, 721, 781, 841]; % SubBlocks are datachunks of 60 trials. Normally, there are 3 subblocks per block. Each subblock will constitute a data point.

elseif BlockSize == 170

    BlockChanges         = [60, 230, 400, 570, 740, 1740 ];                  % (global variables: these variables do not change)
    BlockLengths         = [60, 170, 170, 170, 170, 960  ];
    ProgrammedFirstTrialInBlock    = [1,  61,  231, 401, 571, 741  ];                  % Blocks are session divisions in which  an experimental variable is manipulated (in this case we manipulate RP1 and RP2 blockwise)
    ProgrammedFirstTrialInSubBlock = [1,61,121, 181, 231, 291, 351, 401, 461, 521, 571, 631, 691, 741, 801]; % SubBlocks are datachunks of 60 trials. Normally, there are 3 subblocks per block. Each subblock will constitute a data point.

end

if     Epoch ==1
    RowsOfInterest = find( [Condition(:).LaserDur_ms] == 400 & [Condition(:).Run] == 1 ); % DuringStim Inhibition 400ms
elseif Epoch==2
    RowsOfInterest = find( [Condition(:).LaserDur_ms] == 800 );                           % DuringOut  Inhibition 800ms
end

for iRow=1:length(RowsOfInterest)
    if strcmp(Condition(RowsOfInterest(iRow)).animalID,animalID)
        FinalRow =  RowsOfInterest(iRow);
    end
end

if exist ('FinalRow','var')
else
      error('The selected subject did not run the selected condition') 
end
Schedule =Condition(FinalRow).BiasOrder;          % 1 or 2
Datamatrix = Condition(FinalRow).allBData;



%% Which block has stimulation

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
%% Now lets compute some blockwise metrics of each session

Alldprimes_Sess       = zeros(max(Datamatrix(:,18)),20);
Allcriterions_Sess    = zeros(max(Datamatrix(:,18)),20);
Allcriterions_1_Sess  = zeros(max(Datamatrix(:,18)),20);
AllPR2_Sess           = zeros(max(Datamatrix(:,18)),20);
AllResponseTimeR1     = zeros(max(Datamatrix(:,18)),20);
AllwithdrawalRates    = zeros(max(Datamatrix(:,18)),20);


for iSess = 1:max(Datamatrix(:,18)) % for all sessions

    ThisSessionTrials = find(Datamatrix(:,18)==iSess);
    ThisSessionData   = Datamatrix(ThisSessionTrials,:);
    LastTrial      = length(ThisSessionData);
    ActualFirstTrialInSubBlock =ProgrammedFirstTrialInSubBlock(ProgrammedFirstTrialInSubBlock<LastTrial); % local variable (it depends on the trial number of the session)
    ActualFirstTrialInBlock    =ProgrammedFirstTrialInBlock(ProgrammedFirstTrialInBlock<LastTrial); % local variable (it depends on the trial number of the session)
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

%% Run fisher test to test whether there is am effect of the laser on Pr2
% we will only use the last 60 trials of the to-be-compared blocks (3 and 5)

% SessSplit3    = find(WichBlockIsStimulated==3);
% SessSplit5    = find(WichBlockIsStimulated==5);
% 
% TrialsSplit3_log = logical(sum (Datamatrix(:,18)==SessSplit3,2));
% TrialsSplit5_log = logical(sum (Datamatrix(:,18)==SessSplit5,2));
% 
% MatSplit3 = Datamatrix(TrialsSplit3_log,:);
% MatSplit5 = Datamatrix(TrialsSplit5_log,:);
% 
% [F_Split3]= prepFischertest(MatSplit3,Punishment);
% [h_3,p_3,stats_3] = fishertest(F_Split3);
% 
% [F_Split5]= prepFischertest(MatSplit5,Punishment);
% [h_5,p_5,stats_5] = fishertest(F_Split5);
%% Last variables before Plotting

AllPR2_Sess             = AllPR2_Sess(:,1:length(ActualFirstTrialInSubBlock));
Allcriterions_Sess      = Allcriterions_Sess(:,1:length(ActualFirstTrialInSubBlock));
Allcriterions_1_Sess    = Allcriterions_1_Sess(:,1:length(ActualFirstTrialInSubBlock));
Alldprimes_Sess         = Alldprimes_Sess(:,1:length(ActualFirstTrialInSubBlock));
AllResponseTimeR1       = AllResponseTimeR1(:,1:length(ActualFirstTrialInSubBlock));
AllResponseTimeR2       = AllResponseTimeR2(:,1:length(ActualFirstTrialInSubBlock));
AllWithdrawalRates_Sess = AllWithdrawalRates_Sess(:,1:length(ActualFirstTrialInSubBlock));
AllRewardRates_Sess     = AllRewardRates_Sess(:,1:length(ActualFirstTrialInSubBlock));

if Punishement
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

SubBlockxaxis      = ActualFirstTrialInSubBlock + 29  ;
Blockxaxis         = [30,150,330,510,690,870];
%% Premature withdrawals
LabelConditions = 0;
%% Plotting
orange= [0 0.4470 0.7410];
blue=   [0.8500 0.3250 0.0980];

if Punishement
    Shading_x = [BlockChanges(1),BlockChanges(2);BlockChanges(3),BlockChanges(4)];       %Specify the area of the painted sections. First row block 1, second block 2
    xaxis = [0, 820];
else
    Shading_x = [BlockChanges(2),BlockChanges(3);BlockChanges(4),BlockChanges(5)];       %Specify the area of the painted sections. First row block 1, second block 2
    xaxis = [0, 900];
end
if LabelConditions
    BlockContingenciesString1 = {'.7 : .7','.4 : 1','1 : .4','.4 : 1','1 : .4','.7 : .7'};
    BlockContingenciesString2 = {'.7 : .7','1 : .4','.4 : 1','1 : .4','.4 : 1','.7 : .7'};
    BlockContingenciesString_Pun1 = {'.7 : .7','.7  : .7P','.7 : .7','.7P : .7 ','.7 : .7'};
    BlockContingenciesString_Pun2 = {'.7 : .7','.7P : .7' ,'.7 : .7','.7 : .7 P','.7 : .7'};

    if Punishement
        if Schedule == 1
            BlockContingenciesString = BlockContingenciesString_Pun1 ;
        elseif Schedule ==2
            BlockContingenciesString = BlockContingenciesString_Pun2 ;
        end
    else
        if Schedule == 1
            BlockContingenciesString = BlockContingenciesString1;

        elseif Schedule == 2
            BlockContingenciesString = BlockContingenciesString2;
        end
    end
end

figure

%Shading_x = [240,420;600,780];       %Specify the area of the painted sections. First row block 1, second block 2
basevalue = -2;
if Punishement
    Shading_y = [1,1;3.5,3.5];           %Specify the area of the painted sections. First row block 1, second block 2
    Shading_y_c =[2,2;4,4] ;

else
    Shading_y = [1,1;3,3];           %Specify the area of the painted sections. First row block 1, second block 2
    Shading_y_c =[2,2;4,4] ;
end


ax = arrayfun( @(i) subplot(4,2,i,'NextPlot','add','Box','on'), 1:8 );
%ax = arrayfun( @(i) subplot(2,2,i,'NextPlot','add','Box','on'), [1:4] );
arrayfun( @(i) xline( ax(1), BlockChanges(i)+0.5,'--k'), [1:length(BlockChanges)-1] );                                  % Plot PR2 for Block 3 stimulation
arrayfun( @(i) plot( ax(1), SubBlockxaxis, AllPR2_Sess((WichBlockIsStimulated==OptoInBlocks(1)),:),'color', [.5 .5 .5] ), [1:2] )     % Plot individual sessions
p1 = arrayfun( @(i) plot( ax(1), SubBlockxaxis, MeanPR2(1,:),'-b','LineWidth',3), [1:2] );                                   % Plot mean PR2
arrayfun( @(i) area( ax(1), Shading_x(1,:), Shading_y(1,:),'FaceColor','y','FaceAlpha',.1) , [1:2]);                    % Plot shading stimulated block
if Punishement
    arrayfun( @(i) axis( ax(1), [xaxis 0.15 0.85]), [1:2]);
else
    arrayfun( @(i) axis( ax(1), [xaxis 0.00 1.00]), [1:2]);
end
arrayfun( @(i) xticks(ax(1), 0:200:800), [1:2]);

if LabelConditions
    arrayfun( @(i) text( ax(1),ProgrammedFirstTrialInBlock(1),0.9*ones(1,length(ProgrammedFirstTrialInBlock(1))),BlockContingenciesString(1),'FontSize',14), [1:2]); % adding condition labels in the top. One line for the (1) first, (2) middle and (3) last
    arrayfun( @(i) text( ax(1),ProgrammedFirstTrialInBlock(2:end-1)+60,0.9*ones(1,length(ProgrammedFirstTrialInBlock(2:end-1))),BlockContingenciesString(2:end-1),'FontSize',14), [1:2]);
    arrayfun( @(i) text( ax(1),ProgrammedFirstTrialInBlock(end)+30,0.9*ones(1,length(ProgrammedFirstTrialInBlock(end))),BlockContingenciesString(end),'FontSize',14), 1:2);
end

a = get(gca,'XTickLabel');
b = get(gca,'YTickLabel');
ylabel(ax(1),'PR2','fontsize',12)

arrayfun( @(i) xline( ax(2), BlockChanges(i)+0.5,'--k'), 1:length(BlockChanges)-1 );                                   % Plot PR2 for Block 5 stimulation
arrayfun( @(i) plot( ax(2), SubBlockxaxis, AllPR2_Sess((WichBlockIsStimulated==OptoInBlocks(2)),:),'color', [.5 .5 .5] ), [1:2] );     % Plot individual sessions
arrayfun( @(i) plot( ax(2), SubBlockxaxis, MeanPR2(2,:),'-b' ,'LineWidth',3), [1:2] ) ;                                  % Plot mean PR2
arrayfun( @(i) area( ax(2), Shading_x(2,:), Shading_y(1,:),'FaceColor','y','FaceAlpha',.1) , [1:2]);                     % Plot shading stimulated block
if Punishement
    arrayfun( @(i) axis( ax(2), [xaxis 0.15 0.85]), [1:2]);
else
    arrayfun( @(i) axis( ax(2), [xaxis 0.00 1.00]), [1:2]);
end
if LabelConditions
    arrayfun( @(i) text( ax(2),ProgrammedFirstTrialInBlock(1),0.9*ones(1,length(ProgrammedFirstTrialInBlock(1))),BlockContingenciesString(1),'FontSize',14), [1:2]); % adding condition labels in the top. One line for the (1) first, (2) middle and (3) last
    arrayfun( @(i) text( ax(2),ProgrammedFirstTrialInBlock(2:end-1)+60,0.9*ones(1,length(ProgrammedFirstTrialInBlock(2:end-1))),BlockContingenciesString(2:end-1),'FontSize',14), [1:2]);
    arrayfun( @(i) text( ax(2),ProgrammedFirstTrialInBlock(end)+30,0.9*ones(1,length(ProgrammedFirstTrialInBlock(end))),BlockContingenciesString(end),'FontSize',14), [1:2]);
end

arrayfun( @(i) xline( ax(3), BlockChanges(i)+0.5,'--k'), [1:length(BlockChanges)-1] );                                  % Plot PR2 for Block 3 stimulation
arrayfun( @(i) plot( ax(3), SubBlockxaxis, Allcriterions_Sess((WichBlockIsStimulated==OptoInBlocks(1)),:),'color', [.5 .5 .5] ), [1:2] )     % Plot individual sessions
p1 = arrayfun( @(i) plot( ax(3), SubBlockxaxis, MeanC(1,:),'-b','LineWidth',3), [1:2] );                                   % Plot mean PR2
arrayfun( @(i) area( ax(3), Shading_x(1,:), Shading_y_c(1,:),basevalue,'FaceColor','y','FaceAlpha',.1) , [1:2]);                    % Plot shading stimulated block
if Punishement
    arrayfun( @(i) axis( ax(3), [xaxis -1.1 1.1]), [1:2]);
else
    arrayfun( @(i) axis( ax(3), [xaxis -2 +2]), [1:2]);
end
arrayfun( @(i) xticks(ax(3), 0:200:800), [1:2]);

if LabelConditions
    arrayfun( @(i) text( ax(3),ProgrammedFirstTrialInBlock(1),0.9*ones(1,length(ProgrammedFirstTrialInBlock(1))),BlockContingenciesString(1),'FontSize',14), [1:2]); % adding condition labels in the top. One line for the (1) first, (2) middle and (3) last
    arrayfun( @(i) text( ax(3),ProgrammedFirstTrialInBlock(2:end-1)+60,0.9*ones(1,length(ProgrammedFirstTrialInBlock(2:end-1))),BlockContingenciesString(2:end-1),'FontSize',14), [1:2]);
    arrayfun( @(i) text( ax(3),ProgrammedFirstTrialInBlock(end)+30,0.9*ones(1,length(ProgrammedFirstTrialInBlock(end))),BlockContingenciesString(end),'FontSize',14), [1:2]);
end
%legend (ax(1),[p1(1)], {'R2'},'Location','southwest','FontSize',14)

a = get(gca,'XTickLabel');
b = get(gca,'YTickLabel');
ylabel(ax(3),'c: method 1','fontsize',12)

arrayfun( @(i) xline( ax(4), BlockChanges(i)+0.5,'--k'), [1:length(BlockChanges)-1] );                                   % Plot PR2 for Block 5 stimulation
arrayfun( @(i) plot( ax(4), SubBlockxaxis, Allcriterions_Sess((WichBlockIsStimulated==OptoInBlocks(2)),:),'color', [.5 .5 .5] ), [1:2] );     % Plot individual sessions
arrayfun( @(i) plot( ax(4), SubBlockxaxis, MeanC(2,:),'-b' ,'LineWidth',3), [1:2] ) ;                                  % Plot mean PR2
arrayfun( @(i) area( ax(4), Shading_x(2,:), Shading_y_c(1,:),basevalue,'FaceColor','y','FaceAlpha',.1) , [1:2]);                     % Plot shading stimulated block
if Punishement
    arrayfun( @(i) axis( ax(4), [xaxis -1.1 1.1]), [1:2]);
else
    arrayfun( @(i) axis( ax(4), [xaxis -2 +2]), [1:2]);
end
if LabelConditions
    arrayfun( @(i) text( ax(4),ProgrammedFirstTrialInBlock(1),0.9*ones(1,length(ProgrammedFirstTrialInBlock(1))),BlockContingenciesString(1),'FontSize',14), [1:2]); % adding condition labels in the top. One line for the (1) first, (2) middle and (3) last
    arrayfun( @(i) text( ax(4),ProgrammedFirstTrialInBlock(2:end-1)+60,0.9*ones(1,length(ProgrammedFirstTrialInBlock(2:end-1))),BlockContingenciesString(2:end-1),'FontSize',14), [1:2]);
    arrayfun( @(i) text( ax(4),ProgrammedFirstTrialInBlock(end)+30,0.9*ones(1,length(ProgrammedFirstTrialInBlock(end))),BlockContingenciesString(end),'FontSize',14), [1:2]);
end

arrayfun( @(i) xline( ax(5), BlockChanges(i)+0.5,'--k'), [1:length(BlockChanges)-1] );                                  % Plot PR2 for Block 3 stimulation
arrayfun( @(i) plot( ax(5), SubBlockxaxis, Allcriterions_1_Sess((WichBlockIsStimulated==OptoInBlocks(1)),:),'color', [.5 .5 .5] ), [1:2] )     % Plot individual sessions
p1 = arrayfun( @(i) plot( ax(5), SubBlockxaxis, MeanC_1(1,:),'-b','LineWidth',3), [1:2] );                                   % Plot mean PR2
arrayfun( @(i) area( ax(5), Shading_x(1,:), Shading_y_c(1,:),basevalue,'FaceColor','y','FaceAlpha',.1) , [1:2]);                    % Plot shading stimulated block
if Punishement
    arrayfun( @(i) axis( ax(5), [xaxis -1.1 1.1]), [1:2]);
else
    arrayfun( @(i) axis( ax(5), [xaxis -2 +2]), [1:2]);
end
arrayfun( @(i) xticks(ax(5), 0:200:800), [1:2]);

if LabelConditions
    arrayfun( @(i) text( ax(5),ProgrammedFirstTrialInBlock(1),0.9*ones(1,length(ProgrammedFirstTrialInBlock(1))),BlockContingenciesString(1),'FontSize',14), [1:2]); % adding condition labels in the top. One line for the (1) first, (2) middle and (3) last
    arrayfun( @(i) text( ax(5),ProgrammedFirstTrialInBlock(2:end-1)+60,0.9*ones(1,length(ProgrammedFirstTrialInBlock(2:end-1))),BlockContingenciesString(2:end-1),'FontSize',14), [1:2]);
    arrayfun( @(i) text( ax(5),ProgrammedFirstTrialInBlock(end)+30,0.9*ones(1,length(ProgrammedFirstTrialInBlock(end))),BlockContingenciesString(end),'FontSize',14), [1:2]);
end
%legend (ax(1),[p1(1)], {'R2'},'Location','southwest','FontSize',14)

a = get(gca,'XTickLabel');
b = get(gca,'YTickLabel');
ylabel(ax(5),'c: method 2','fontsize',12)

arrayfun( @(i) xline( ax(6), BlockChanges(i)+0.5,'--k'), [1:length(BlockChanges)-1] );                                   % Plot PR2 for Block 5 stimulation
arrayfun( @(i) plot( ax(6), SubBlockxaxis, Allcriterions_1_Sess((WichBlockIsStimulated==OptoInBlocks(2)),:),'color', [.5 .5 .5] ), [1:2] );     % Plot individual sessions
arrayfun( @(i) plot( ax(6), SubBlockxaxis, MeanC_1(2,:),'-b' ,'LineWidth',3), [1:2] ) ;                                  % Plot mean PR2
arrayfun( @(i) area( ax(6), Shading_x(2,:), Shading_y_c(1,:),basevalue,'FaceColor','y','FaceAlpha',.1) , [1:2]);                     % Plot shading stimulated block
if Punishement
    arrayfun( @(i) axis( ax(6), [xaxis -1.1 1.1]), [1:2]);
else
    arrayfun( @(i) axis( ax(6), [xaxis -2 2]), [1:2]);
end
if LabelConditions
    arrayfun( @(i) text( ax(6),ProgrammedFirstTrialInBlock(1),0.9*ones(1,length(ProgrammedFirstTrialInBlock(1))),BlockContingenciesString(1),'FontSize',14), [1:2]); % adding condition labels in the top. One line for the (1) first, (2) middle and (3) last
    arrayfun( @(i) text( ax(6),ProgrammedFirstTrialInBlock(2:end-1)+60,0.9*ones(1,length(ProgrammedFirstTrialInBlock(2:end-1))),BlockContingenciesString(2:end-1),'FontSize',14), [1:2]);
    arrayfun( @(i) text( ax(6),ProgrammedFirstTrialInBlock(end)+30,0.9*ones(1,length(ProgrammedFirstTrialInBlock(end))),BlockContingenciesString(end),'FontSize',14), [1:2]);
end


arrayfun( @(i) xline( ax(7), BlockChanges(i)+0.5,'--k'), [1:length(BlockChanges)-1] );                                   % Plot dprimes for Block 3 stimulation
%arrayfun( @(i) plot( ax(7), SubBlockxaxis, Alldprimes_Sess((WichBlockIsStimulated==OptoInBlocks(1)),:),'color', [.5 .5 .5] ), [1:2] ); % Plot individual sessions
p3 = arrayfun( @(i) plot( ax(7), SubBlockxaxis, MeanRewardRates(1,:),'-k','LineWidth',3), [1:2] );
p4 = arrayfun( @(i) plot( ax(7), SubBlockxaxis, MeanWithdrawalRates(1,:),'-g','LineWidth',3), [1:2] );
arrayfun( @(i) area( ax(7), Shading_x(1,:), Shading_y(2,:),'FaceColor','y','FaceAlpha',.1) , [1:2]);                     % Plot shading stimulated block
arrayfun( @(i) axis (ax(7),[xaxis 0 1]), [1:2]);

ylabel(ax(7),'Prob.')
xlabel(ax(7),'Trial Nr')

legend (ax(7),[p3(1) p4(1)], 'Reward','Abort','Location','southwest','FontSize',14)


arrayfun( @(i) xline( ax(8), BlockChanges(i)+0.5,'--k'), [1:length(BlockChanges)-1] );                                   % Plot dprimes for Block 5 stimulation
arrayfun( @(i) plot( ax(8), SubBlockxaxis, MeanRewardRates(2,:),'-k' ,'LineWidth',3), [1:2] );
arrayfun( @(i) plot( ax(8), SubBlockxaxis, MeanWithdrawalRates(2,:),'-g','LineWidth',3), [1:2] );
arrayfun( @(i) area( ax(8), Shading_x(2,:), Shading_y(2,:),'FaceColor','y','FaceAlpha',.1) , [1:2]);                     % Plot shading stimulated block
arrayfun( @(i) axis (ax(8),[xaxis 0 1]), [1:2]);

xlabel(ax(8),'Trial Nr')





%% All opto inactivation)

if ~ Punishement
Block_2_datapoints = 1:6;
Block_3_datapoints = 4:8;
Block_4_datapoints = 6:12;
Block_5_datapoints = 10:14;

elseif Punishement
Block_2_datapoints = 2:7;
Block_3_datapoints = 6:9;
Block_4_datapoints = 8:13;
PR2_Block2_ON  = PR2_Block3_ON;
PR2_Block4_ON  = PR2_Block5_ON;
c_Block2_ON    = c_Block3_ON;
c_Block4_ON    = c_Block5_ON;
clearvars PR2_Block3_ON PR2_Block5_ON c_Block3_ON c_Block5_ON 
end


if ~ Punishement
PR2_Laser_ON_All = [PR2_Block3_ON(:,Block_3_datapoints);PR2_Block5_ON(:,Block_5_datapoints)];
PR2_Laser_ON_All(end+1,:)=mean(PR2_Laser_ON_All,1);
PR2_Laser_OFF_All = [PR2_Block3_ON(:,Block_5_datapoints);PR2_Block5_ON(:,Block_3_datapoints)];
PR2_Laser_OFF_All(end+1,:)=mean(PR2_Laser_OFF_All,1);

c_Laser_ON_All = [c_Block3_ON(:,Block_3_datapoints);c_Block5_ON(:,Block_5_datapoints)];
c_Laser_ON_All(end+1,:)=mean(c_Laser_ON_All,1);
c_Laser_OFF_All = [c_Block3_ON(:,Block_5_datapoints);c_Block5_ON(:,Block_3_datapoints)];
c_Laser_OFF_All(end+1,:)=mean(c_Laser_OFF_All,1);
else

PR2_Laser_ON_All = [PR2_Block2_ON(:,Block_2_datapoints);PR2_Block4_ON(:,Block_4_datapoints)];
PR2_Laser_OFF_All = [PR2_Block2_ON(:,Block_4_datapoints);PR2_Block4_ON(:,Block_2_datapoints)];
c_Laser_ON_All = [c_Block2_ON(:,Block_2_datapoints);c_Block4_ON(:,Block_4_datapoints)];
c_Laser_OFF_All = [c_Block2_ON(:,Block_4_datapoints);c_Block4_ON(:,Block_2_datapoints)];
end
 
Block_relevant_changes_B2 =[BlockChanges(1:2)];
Block_relevant_changes_B3 =[BlockChanges(2:3)];
Block_relevant_changes_B4 =[BlockChanges(3:4)];
Block_relevant_changes_B5 =[BlockChanges(4:5)]; 
MeanPR2_Laser_ON_Block2  = MeanPR2(1,Block_2_datapoints);
MeanPR2_Laser_OFF_Block2 = MeanPR2(2,Block_2_datapoints);

figure
if ~Punishement
Blocks_to_plot = Block_3_datapoints;
Trialzero = 240;
Block_edges = Block_relevant_changes_B3-Trialzero;
elseif Punishement
Blocks_to_plot = Block_2_datapoints;
Trialzero = 120;
Block_edges = Block_relevant_changes_B2-Trialzero;
end

ax = arrayfun( @(i) subplot(1,1,i,'NextPlot','add','Box','on'), [1:1] );
    arrayfun( @(i) xline( ax(1), Block_edges(i)+0.5,'--k'), [1:length(Block_edges)] );                                  % Plot PR2 for Block 3 stimulation

if Punishement
    arrayfun( @(i) plot( ax(1), SubBlockxaxis(Blocks_to_plot)-Trialzero, c_Laser_ON_All(1:end,:),"--",'color', [0.9294,0.6941,0.1255],'LineWidth',0.5 ), [1] )     % Plot individual sessions LASER ON
    arrayfun( @(i) plot( ax(1), SubBlockxaxis(Blocks_to_plot)-Trialzero, c_Laser_OFF_All(1:end,:),"--",'color', [0,0,0],'LineWidth',0.01 ), [1] )     % Plot individual sessions LASER OFF
    p1 = arrayfun( @(i) plot( ax(1), SubBlockxaxis(Blocks_to_plot)-Trialzero, mean (Allcriterions_1_Sess(WichBlockIsStimulated==2,Block_2_datapoints),1),'color',[0.9294,0.6941,0.1255],'LineWidth',3), [1] );                                              % Plot mean PR2 LASER ON
         arrayfun( @(i) plot( ax(1), SubBlockxaxis(Blocks_to_plot)-Trialzero, mean (Allcriterions_1_Sess(WichBlockIsStimulated==4,Block_4_datapoints),1),'color',[0.9294,0.6941,0.1255],'LineWidth',3), [1] );                                              % Plot mean PR2 LASER ON
    p2 = arrayfun( @(i) plot( ax(1), SubBlockxaxis(Blocks_to_plot)-Trialzero,mean (Allcriterions_1_Sess(WichBlockIsStimulated==4,Block_2_datapoints),1),'color',[0,0,0],'LineWidth',3), [1] );                                              % Plot mean PR2 LASER OFF
         arrayfun( @(i) plot( ax(1), SubBlockxaxis(Blocks_to_plot)-Trialzero,mean (Allcriterions_1_Sess(WichBlockIsStimulated==2,Block_4_datapoints),1),'color',[0,0,0],'LineWidth',3), [1] );                                              % Plot mean PR2 LASER OFF
    
    arrayfun( @(i) axis( ax(1), [[SubBlockxaxis(Blocks_to_plot(1))-Trialzero,SubBlockxaxis(Blocks_to_plot(end))-Trialzero] -1.1 1.1]), [1]);

else
    arrayfun( @(i) plot( ax(1), SubBlockxaxis(Blocks_to_plot)-Trialzero, c_Laser_ON_All(1:end,:),"--",'color', [0.9294,0.6941,0.1255],'LineWidth',0.5 ), [1] )     % Plot individual sessions LASER ON
    arrayfun( @(i) plot( ax(1), SubBlockxaxis(Blocks_to_plot)-Trialzero, c_Laser_OFF_All(1:end,:),"--",'color', [0,0,0],'LineWidth',0.01 ), [1] )     % Plot individual sessions LASER OFF
    p1 = arrayfun( @(i) plot( ax(1), SubBlockxaxis(Blocks_to_plot)-Trialzero, c_Laser_ON_All(end,:),'color',[0.9294,0.6941,0.1255],'LineWidth',3), [1] );                                              % Plot mean PR2 LASER ON
    p2 = arrayfun( @(i) plot( ax(1), SubBlockxaxis(Blocks_to_plot)-Trialzero,c_Laser_OFF_All(end,:),'color',[0,0,0],'LineWidth',3), [1] );                                              % Plot mean PR2 LASER OFF

    arrayfun( @(i) axis( ax(1), [[SubBlockxaxis(Blocks_to_plot(1))-Trialzero,SubBlockxaxis(Blocks_to_plot(end))-Trialzero] -2.000 2.00]), [1]);
end
arrayfun( @(i) xticks(ax(1), 0:200:800), [1]);

if LabelConditions
    arrayfun( @(i) text( ax(1),ProgrammedFirstTrialInBlock(1),0.9*ones(1,length(ProgrammedFirstTrialInBlock(1))),BlockContingenciesString(1),'FontSize',14), [1:2]); % adding condition labels in the top. One line for the (1) first, (2) middle and (3) last
    arrayfun( @(i) text( ax(1),ProgrammedFirstTrialInBlock(2:end-1)+60,0.9*ones(1,length(ProgrammedFirstTrialInBlock(2:end-1))),BlockContingenciesString(2:end-1),'FontSize',14), [1:2]);
    arrayfun( @(i) text( ax(1),ProgrammedFirstTrialInBlock(end)+30,0.9*ones(1,length(ProgrammedFirstTrialInBlock(end))),BlockContingenciesString(end),'FontSize',14), [1:2]);
end
a = get(gca,'XTickLabel');
b = get(gca,'YTickLabel');
ylabel(ax(1),'c','fontsize',12)
title (ax(1),'Lights ON vs Lights OFF')
axis square
axis (ax(1),'square')
xlabel('Trials from block change')
legend (ax(1),[p1 p2], 'ON','OFF','Location','southwest','FontSize',14)

%% HELPER FUNCTIONS

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


