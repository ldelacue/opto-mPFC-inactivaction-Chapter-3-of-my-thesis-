%% GLM_Fit_real_data.m

clear all
load("SPP_Blockwise_opto_struct.mat")  %alternatively you can upload SPP 
Condition = SPP_Blockwise_opto_struct; %alternatively you can upload SPP 
Epoch     = 1; % 1 During Stim 2 Outcome (no outcome dataset for SPP manipualtions)
Dataframe = 0; % Should be 0 for historical reasons, but has no importance here.
TrialsInPast = [1,2,3,4];
RowsOfInterest_Halo = [];
RowsOfInterest_YFP =  [];

if     Epoch ==1
RowsOfInterest = find( [Condition(:).LaserDur_ms] == 400 & [Condition(:).Run] == 1 ); % DuringStim Inhibition 400ms
elseif Epoch==2                                                                          
RowsOfInterest = find( [Condition(:).LaserDur_ms] == 800 );                           % DuringOut  Inhibition 800ms
end

for iRow=1:length(RowsOfInterest)
    if strcmp(Condition(RowsOfInterest(iRow)).Pheno,'Halo')
      RowsOfInterest_Halo =  [RowsOfInterest_Halo, RowsOfInterest(iRow)];
    else strcmp(Condition(RowsOfInterest(iRow)).Pheno,'YFP')
      RowsOfInterest_YFP =  [RowsOfInterest_YFP, RowsOfInterest(iRow)];
    end
end
RowsOfInterest = RowsOfInterest_Halo;                       % CHOOSE PHENOTYPE (RowsOfInterest_Halo or RowsOfInterest_YFP)
for iDataset=1:length(RowsOfInterest)

    DataMatrixforIndexes = Condition(RowsOfInterest(iDataset)).allBData ;
    DataMatrix = [Condition(RowsOfInterest(iDataset)).allBDataKDB,DataMatrixforIndexes(:,end)] ;
    DataMatrixforIndexes = DataMatrixforIndexes(~isnan(DataMatrixforIndexes(:,5)),:);
    DataMatrix = DataMatrix(~isnan(DataMatrix(:,3)),:);
    BlocksofInterest = [3,5];

    NewDataMatrixforIndexes    = DataMatrixforIndexes(DataMatrixforIndexes(:,17)==BlocksofInterest(1)|DataMatrixforIndexes(:,17)==BlocksofInterest(2),:); % Select blocks of interest
    NewDataMatrix         = DataMatrix(DataMatrixforIndexes(:,17)==BlocksofInterest(1)|DataMatrixforIndexes(:,17)==BlocksofInterest(2),:);
    LaserONDataMatrix     = NewDataMatrix(NewDataMatrixforIndexes(:,15)==1,:);
    LaserOFFDataMatrix    = NewDataMatrix(NewDataMatrixforIndexes(:,15)==0,:);
    LaserON_Ind           = find(NewDataMatrixforIndexes(:,15)==1);
    LaserOFF_Ind          = find(NewDataMatrixforIndexes(:,15)==0);
    OneHotEncoder    = 1;  % 1 for Session by Session Regressor and 2 for Condition by Condition Regressor
    [Regressors,DependentVar, betas, pvalues, SEWeights, VIF] = fitGLM_opto_on_off(NewDataMatrix,LaserON_Ind,LaserOFF_Ind,OneHotEncoder,Dataframe,TrialsInPast);

    FittedStimulusMeans(iDataset,:)           = betas(1:6);
    RewardedWeights_LON(iDataset,:)           = betas(7:(7+length(TrialsInPast)-1) );  %After rewarded trial
    RewardedWeights_LOFF(iDataset,:)          = betas(7+length(TrialsInPast) : (7+ 2*length(TrialsInPast)-1) );  %After rewarded trial
    UnrewardedWeights_LON(iDataset,:)         = betas(7+2*length(TrialsInPast) : (7+ 3*length(TrialsInPast)-1) ); %After unrewarded trial
    UnrewardedWeights_LOFF(iDataset,:)         = betas(7+3*length(TrialsInPast) : (7+ 4*length(TrialsInPast)-1) ); %After unrewarded trial

end

    RewardedWeights_LON(size(RewardedWeights_LON,1)+1,:)       = mean(RewardedWeights_LON,1);
    UnrewardedWeights_LON(size(UnrewardedWeights_LON,1)+1,:)   = mean(UnrewardedWeights_LON,1);

    RewardedWeights_LOFF(size(RewardedWeights_LOFF,1)+1,:)     = mean(RewardedWeights_LOFF,1);
    UnrewardedWeights_LOFF(size(UnrewardedWeights_LOFF,1)+1,:) = mean(UnrewardedWeights_LOFF,1);

    %% PLotting

figure 
ax1 = subplot(1,3,1);
xaxis = TrialsInPast;
    for iDataset=1:length(RowsOfInterest)
        p1 = plot(xaxis(~isnan(RewardedWeights_LON(iDataset,:))),RewardedWeights_LON(iDataset,~isnan(RewardedWeights_LON(iDataset,:))),'b--');
        hold on;
        plot(xaxis(~isnan(RewardedWeights_LON(end,:))),RewardedWeights_LON(end,~isnan(RewardedWeights_LON(end,:))),'b-','LineWidth',2)
        hold on;

        p2 = plot(xaxis(~isnan(UnrewardedWeights_LON(iDataset,:))),UnrewardedWeights_LON(iDataset,~isnan(UnrewardedWeights_LON(iDataset,:))),'r--');
        hold on;
        plot(xaxis(~isnan(UnrewardedWeights_LON(end,:))),UnrewardedWeights_LON(end,~isnan(UnrewardedWeights_LON(end,:))),'r-','LineWidth',2)
        hold on
    end

    xlabel('Trials in the past')
    ylabel('Weights')
    title(ax1,'Influence of outcome on current choice - Laser ON trials')
    legend([p1 p2],{'Reward', 'TimeOut'})
    axis square
    ylim ([-0.2, 1.0])

ax2 = subplot(1,3,2);

    for iDataset=1:length(RowsOfInterest)
        p1 = plot(xaxis(~isnan(RewardedWeights_LOFF(iDataset,:))),RewardedWeights_LOFF(iDataset,~isnan(RewardedWeights_LON(iDataset,:))),'b--');
        hold on;
        plot(xaxis(~isnan(RewardedWeights_LOFF(end,:))),RewardedWeights_LOFF(end,~isnan(RewardedWeights_LON(end,:))),'b-','LineWidth',2)
        hold on;

        p2 = plot(xaxis(~isnan(UnrewardedWeights_LOFF(iDataset,:))),UnrewardedWeights_LOFF(iDataset,~isnan(UnrewardedWeights_LOFF(iDataset,:))),'r--');
        hold on;
        plot(xaxis(~isnan(UnrewardedWeights_LOFF(end,:))),UnrewardedWeights_LOFF(end,~isnan(UnrewardedWeights_LOFF(end,:))),'r-','LineWidth',2)
        hold on
    end
    xlabel('Trials in the past')
    ylabel('Weights')
    title(ax2,'Influence of outcome history on current choice - Laser OFF trials')
    legend([p1 p2],{'Reward', 'TimeOut'})
    axis square
    ylim([-0.2, 1.0])

    ax3 = subplot(1,3,3);

legend([p1 p2],{'Reward', 'TimeOut'})
axis square
ylim([-0.2, 1.0])
xlim([0.75, 4.25])

subplot(1,3,3);
for iDataset=1:length(RowsOfInterest)
    p1 = plot(xaxis(~isnan(RewardedWeights_LON(iDataset,:))),RewardedWeights_LON(iDataset,~isnan(RewardedWeights_LON(iDataset,:))),'--','Color',[0.9290 0.6940 0.1250]);
    hold on;
    plot(xaxis(~isnan(RewardedWeights_LON(end,:))),RewardedWeights_LON(end,~isnan(RewardedWeights_LON(end,:))),'-','LineWidth',2,'Color',[0.9290 0.6940 0.1250])
    hold on;

    p2 = plot(xaxis(~isnan(RewardedWeights_LOFF(iDataset,:))),RewardedWeights_LOFF(iDataset,~isnan(RewardedWeights_LOFF(iDataset,:))),'k--');
    hold on;
    plot(xaxis(~isnan(RewardedWeights_LOFF(end,:))),RewardedWeights_LOFF(end,~isnan(RewardedWeights_LOFF(end,:))),'k-','LineWidth',2)
    hold on;
end
legend([p1 p2],{'Laser ON', 'Laser OFF'})
axis square
ylim([-0.2, 1.0])
xlim([ 0.75, 4.25])
title('Influence of reward history on current choice')


 
