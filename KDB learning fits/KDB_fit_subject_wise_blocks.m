%% KDB_fit_subject_wise_blocks

% This code reproduces figures 3.8 and S3.4 of my thesis. It plots
% parameters returned by the kdb learning models for either RP or SPP
% experiments in Laser ON / OFF blocks during stimulus (Epoch=1) or during
% outcome (Epoch = 2).
% 
% Needs kdbfit3x.m and BIC_RP_Stim, BIC_RP_Out, BIC_RP_Stim_Run2,
% BIC_SPP_Stim, KDB_Fitted_params_Project1.m.
% 
% Note: In my local kdbfit3x.m version, I restricted the range of the gamma parameter to
% 0.7-1. This helps the gamma parameter to be constrained
% to a reasonable range, in instances where there are little trials in the
% steady state (as it is the case in these datasets, where contingencies are
% manipualted several times per session (on average here: every 150 completed
% trials)). I also set plotting = 0 within kdbfit3x.m.
% You would have to the same in your kdbfit3x.m version 
% if you want to reproduce the exact same results as in the Chapter 3 of my thesis.
% 
% I used model IR-RD (in de la Cuesta et al., 2025), but with an additional learning rate
% for blocks in which the laser was ON. This is implemented in kdbfit3x.m by inputing "model 3b" and
% specifying the delta_mapping variable blockwise. So do not be confused, in spirit, we
% are still fitting model 3a, one learning rate for all stimuli, only that
% since we have an additional learning rate for Laser ON trials, model 3b
% instead of model 3a is the input to kdbfit3x.m (see documentation kdbfit3x.m). 

%                                               Luis de la Cuesta 02/2025

clear all  
close all

Manipulation = 'RP';
Epoch     = 1; % 1 During Stim 2 Outcome (only in RP experiment)
Run_nr    = 1;
%load('StimulusMeans.mat')

if strcmp(Manipulation,'RP')
    load ("RP_Blockwise_opto_struct.mat")
    Condition = RP_Blockwise_opto_struct;
    if Epoch ==1
    ylim = [0.05 0.40];
    elseif Epoch ==2
    ylim = [0.1 0.4];
    end
elseif strcmp(Manipulation,'SPP')
    load ("SPP_Blockwise_opto_struct.mat")
    Condition = SPP_Blockwise_opto_struct;
    ylim = [0.1 0.45];
end

if     Epoch ==1
    RowsOfInterest = find( [Condition(:).LaserDur_ms] == 400 & [Condition(:).Run] == Run_nr); % DuringStim Inhibition 400ms
elseif Epoch==2
    RowsOfInterest = find( [Condition(:).LaserDur_ms] == 800 );                               % DuringOut  Inhibition 800ms
end
Pheno = zeros(1,length(RowsOfInterest));
for iDataset=1:length(RowsOfInterest)
      iRow                   = RowsOfInterest(iDataset);
      Datamatrix_KDB = Condition(iRow).allBDataKDB;
      Datamatrix     = Condition(iRow).allBData;
   %   Datamatrix = Datamatrix(Datamatrix(:,18)==1,:);                % delete. This is only for getting single session fit plots
   %   Datamatrix_KDB = Datamatrix_KDB(Datamatrix(:,18)==1,:);        % delete. This is only for getting single session fit plots
      if strcmp(Condition(iRow).Pheno,'Halo');Pheno(iDataset)=1;end
      Opto_column = Datamatrix(:,15);

    % Which block has stimulation

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

Cum_Block_column = Datamatrix(:,17); % Will be provided later on as Condition
Cum_SubBlock_column   = Datamatrix(:,18); % Will be provided later on as Session

Block_opto_aggr = zeros(max(Datamatrix(:,end)),6);
for iSess= 1:size(Block_opto_aggr,1)
Block_opto_aggr(iSess,WichBlockIsStimulated(iSess))=1;
if iSess>1
    ThisSessionTrialInd = Datamatrix( Datamatrix(:,18)==iSess,1);
    LastSessionBlockInd = Cum_Block_column( Datamatrix(:,18)==iSess-1);
    ThisSessionBlockInd = Cum_Block_column( Datamatrix(:,18)==iSess);
    Cum_Block_column(Datamatrix(:,18)==iSess) = ThisSessionBlockInd+LastSessionBlockInd(end,1);
end
end
LearningRatesCode = Block_opto_aggr+1;
%LearningRatesCode = ones(max(Datamatrix(:,end)),6); % for one learning rate instead of 2
LearningRatesCode = LearningRatesCode';
LearningRatesCode = LearningRatesCode(:);
ConditionCodes    = (1:1:size(Block_opto_aggr,1)*6)';

Datamatrix_KDB = [Datamatrix_KDB(:,1:5),Cum_Block_column,Cum_Block_column];
trialgroup = 'condition';
delta_mapping = dictionary(ConditionCodes', LearningRatesCode'); 
no_update_no_pullback =1;
%no_update_no_pullback =0;
Model      = 'model 3b';
[w_3a_opto(iDataset,:),nll,A,param_3a_opto(iDataset,:),cval,pS2,BIC_3a_opto_RD(iDataset)]= kdbfit3x(Datamatrix_KDB,Model,[],[],[],[],[],'shuffled',no_update_no_pullback,delta_mapping,trialgroup);
clear Cum_Block_column Block_opto_aggr
end
close all
%% Plotting
figure
w = w_3a_opto;
param = param_3a_opto;
orange= [1 0.7 0.3];
ax1 = subplot(1,3,1);
title('Halo')
boxplot(w(Pheno==1,7:8),'Labels',{'delta','delta opto'});hold on
scatter(ones(length(w(Pheno==1,7)),1),w(Pheno==1,7),'x');hold on;scatter(2*ones(length(w(Pheno==1,8)),1),w(Pheno==1,8),'x')
axis square
axis([0 3 ylim])
title(ax1,'HaloR')
ax2 = subplot(1,3,2);
title('YFP')
boxplot(w(Pheno==0,7:8),'Labels',{'delta','delta opto'});hold on
scatter(ones(length(w(Pheno==0,7)),1),w(Pheno==0,7),'x');hold on;scatter(2*ones(length(w(Pheno==0,8)),1),w(Pheno==0,8),'x')
axis square
axis([0 3 ylim])
title(ax2,'YFP')

subplot(1,3,3)
x_OFF = 1-param(Pheno==1);
x2 = 1-param(Pheno==0);
x  = [x_OFF;x2];
% g1 = repmat({'gamma Halo'},size(w(Pheno==1,7:8),1),1);
% g2 = repmat({'gamma YFP'}, size(w(Pheno==0,7:8),1),1);

g1 = repmat({'gamma Halo'},length(param(Pheno==1)),1 );
g2 = repmat({'gamma YFP'}, length(param(Pheno==0)),1 );
g  = [g1; g2];
boxplot(x,g);hold on
scatter(ones(length(param(Pheno==1)),1),1-param(Pheno==1),'x');hold on;scatter(2*ones(length(param(Pheno==0)),1),1-param(Pheno==0),'x')
axis square
axis([0 3 ylim])

% stats

[h,p] = ttest(w(Pheno==0,7), w(Pheno==0,8));

%% Compare BIC between models with and without additional learning rates in Halo vs YFPs subjects


if strcmp(Manipulation,'RP') && Run_nr==1 && Epoch==1
    load('BICs_RP_Stim.mat')

elseif strcmp(Manipulation,'RP') && Run_nr==1 && Epoch==2
    load('BICs_RP_Out.mat')

elseif strcmp(Manipulation,'RP') && Run_nr==2 && Epoch==1
    load('BICs_RP_Stim_Run2.mat')

elseif strcmp(Manipulation,'SPP') && Run_nr==1 && Epoch==1
    load ("SPP_Blockwise_opto_struct.mat")

end
Labels = {'Integrate Rewards (3a)','Integrate Rewards - Opto (3a)'};
figure
ax1 = subplot(1,2,1);
BICs = [BIC_3a_RD- BIC_3a_opto_RD; BIC_3a_opto_RD - BIC_3a_opto_RD];
boxplot(BICs(:,(Pheno==1))','Labels',Labels); hold on
plot( BICs(:,(Pheno==1)),'k-x');hold on
line([-10,10],[0 0],'LineStyle','--')
axis([ 0.5 2.5 -20 20])
axis square
title(ax1,'HaloR')

ax2 = subplot(1,2,2);
boxplot(BICs(:,(Pheno==0))','Labels',Labels); hold on
plot( BICs(:,(Pheno==0)),'k-x');hold on
line([-10,10],[0 0],'LineStyle','--')
axis([ 0.5 2.5 -20 20])
axis square
title(ax2,'YFP')



%% Plot gamma vs delta correlation in all fitted datasets in this project (opto) and of the subjects in de la Cuesta 2025., with model IR-RD
% generate Figure S3.4 chapter 3 of my thesis
load('KDB_Fitted_params_Project1.mat')
figure
w = w_3a_opto;
param = param_3a_opto;
p1 = scatter(1-param(Pheno==1),w(Pheno==1,7),70,'k','x'); hold on
p2 = scatter(1-param(Pheno==1),w(Pheno==1,8),70, orange,'x'); hold on
p3 = scatter(1-KDB_Fitted_Parameters_Long_3a_RD(:,end),KDB_Fitted_Parameters_Long_3a_RD(:,end-1),50,'green','x');hold on
p4 = scatter(1-KDB_Fitted_Parameters_Pigeons_3a_RD(:,end),KDB_Fitted_Parameters_Pigeons_3a_RD(:,end-1),50,'green','o');
ylabel('delta')
xlabel('1-gamma')
axis square
axis([0 0.40 0 0.40])
legend ([p1,p2,p3,p4],{'HaloR - Laser ON','HaloR - Laser OFF','Rats de la Cuesta et al., 2025','Pigeons de la Cuesta et al., 2025'})


