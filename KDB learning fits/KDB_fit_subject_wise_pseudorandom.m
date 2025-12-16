%% KDB_fit_subject_wise_pseudorandom

% This code reproduces figures S3.4D,E,F of my thesis. It plots
% parameters returned by the kdb learning models Laser ON / OFF trials (during stim).

% Needs kdbfit3x_beta.m and BIC_Pseudorandom_Stim.

% Note: kdbfit3x_beta.m is a new version of kdbfit3x.m in which I 
% implemented new functionality to fit different learning rates in Laser ON vs OFF trials 
% that happen in alternation such us in the
% pseudorandom Laser presentation schedule. So far the code could not do
% this. I added this functionality by adding an additional input that specifies how
% trials and learning rates are matched, i.e., trialgroup_def. See kdbfit3x_beta.m

%                                               Luis de la Cuesta 02/2025


clear all  
close all
Exclude_first_subjects     = 1;
% We excluded the first subjects because they underwent testing on the pseudorandom condition 
% after having gone through the blocked RP manipulations sessions
% and we found that animals seemed to have associated optic stimulation
% with the response associated to the contingency with which it was
% previously associated in the previous experiment (RP condition). We decided to not interpret
% the data of these initial animals for these reason. In the rest of
% animals, they were allways tested before in the pseudorandom condition
% prior to the blocked RP manipulations. 

    load ("Pseudorandom_opto_struct.mat")
    Condition = PseudorandomOptoStruct_beta;

RowsOfInterest = 1:length(Condition); % DuringStim Inhibition 400ms
if     Exclude_first_subjects ==1
    RowsOfInterest = RowsOfInterest(4:end);
end

Pheno = zeros(1,length(RowsOfInterest));
for iDataset=1:length(RowsOfInterest)
      iRow                   = RowsOfInterest(iDataset);
      Datamatrix_KDB = Condition(iRow).allBDataKDB;
      Datamatrix     = Condition(iRow).allBData;
      if strcmp(Condition(iRow).Pheno,'Halo');Pheno(iDataset)=1;end
      Opto_column = Datamatrix(:,15);

    % Which block has stimulation

%     for iSess = 1:max(Datamatrix(:,18))
%         SessionTrialInd              = find( Datamatrix(:,18) == iSess);
%         SessionOptoTrialInd          = find( Datamatrix(SessionTrialInd,15)==1);
%         if isempty(SessionOptoTrialInd) % There was not any stimulation in any block!
%             Mssg =  "There was no laser stimulation in the session Nr"+ iSess;
%             error(Mssg)
%         end
%         WichBlockIsStimulated(iSess) = Datamatrix(SessionOptoTrialInd(1),17);
%     end
    clear iSess

Datamatrix_KDB = [Datamatrix_KDB,Datamatrix_KDB(:,6)];
LearningRatesCode = Opto_column+1;
%LearningRatesCode = ones(max(Datamatrix(:,end)),6); % for one learning rate instead of 2
LearningRatesCode = LearningRatesCode';
LearningRatesCode = LearningRatesCode(:);

trialgroup     = 'undefined'; 
trialgroup_def = LearningRatesCode;
delta_mapping = dictionary([1 2], [1 2]);
no_update_no_pullback =1;
Model      = 'model 3b';
[w_3a_opto(iDataset,:),nll,A,param_3a_opto(iDataset,:),cval,pS2,BIC_3a_opto_RD(iDataset)]= kdbfit3x_beta(Datamatrix_KDB,Model,[],[],[],[],[],'shuffled',no_update_no_pullback,delta_mapping,trialgroup,trialgroup_def);
close all
end
%%
figure
w = w_3a_opto;
param = param_3a_opto;
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


subplot(1,3,3) % gammas look very weird.
x_OFF = 1-param(Pheno==1);
x2 = 1-param(Pheno==0);
x  = [x_OFF;x2];
g1 = repmat({'gamma Halo'},length(param(Pheno==1)),1 );
g2 = repmat({'gamma YFP'}, length(param(Pheno==0)),1 );
g  = [g1; g2];
boxplot(x,g);hold on
scatter(ones(length(param(Pheno==1)),1),1-param(Pheno==1),'x');hold on;scatter(2*ones(length(param(Pheno==0)),1),1-param(Pheno==0),'x')
axis square
axis([0 3 ylim])

% stats

[h,p] = ttest(w(Pheno==1,7), w(Pheno==1,8));

%%

load('BICs_Pseudorandom_Stim.mat')
Labels = {'Integrate Rewards (3a)','Integrate Rewards - Opto (3a)'};
figure
subplot(1,2,1)
BICs = [BIC_3a_RD- BIC_3a_opto_RD; BIC_3a_opto_RD - BIC_3a_opto_RD];
boxplot(BICs(:,(Pheno==1))','Labels',Labels); hold on
plot( BICs(:,(Pheno==1)),'k-x');hold on
line([-10,10],[0 0],'LineStyle','--')
axis([ 0.5 2.5 -20 20])
axis square

subplot(1,2,2)
boxplot(BICs(:,(Pheno==0))','Labels',Labels); hold on
plot( BICs(:,(Pheno==0)),'k-x');hold on
line([-10,10],[0 0],'LineStyle','--')
axis([ 0.5 2.5 -20 20])
axis square




