%% GLM_Simulate_hypotheses

% Simulate different hypotheses with different kdb parameters (model IR-RD)
% and fit a GLM model to assess how reward are integrated in each of the
% different algorithmic hypotheses to provide a benchmark to the real data,
% and help decide wich of the hypotheses is the most likely. This generates
% Figure 3.9 of my thesis (simulated data only), for real data go to code
% (GLM_Fit_real_data.m). SImulations happen through msimCDExperiment.m
% code. I implemented some small changes. Just in case, use the version I
% upload.

%                                                   Luis de la Cuesta
%                                                   02/2025

clear all

% For GLM fit
TrialsInPast = [1,2,3,4];%,5,6,7,8];
OneHotEncoder = 1; % Session-wise regressors

% For simulations
nr_sim      = 8;                          % as the number of Halo subjects in the experiment
in.d        = [-1.5,-1,-0.5,0.5,1,1.5];   % stimulus means
in.category = [1,1,1,2,2,2];    % stimulus category
in.model    = 'income RD';         % try 'income', 'unrewarded', 'error', and 'general', 'Lak', 'drift', and 'fixed'
in.preinf   = [1,1,1,0.4,0.4,0.4; 0.4,0.4,0.4,1,1,1; 1,1,1,0.4,0.4,0.4; 0.4,0.4,0.4,1,1,1; 1,1,1,0.4,0.4,0.4; 0.4,0.4,0.4,1,1,1;...
               1,1,1,0.4,0.4,0.4; 0.4,0.4,0.4,1,1,1; 1,1,1,0.4,0.4,0.4; 0.4,0.4,0.4,1,1,1; 1,1,1,0.4,0.4,0.4; 0.4,0.4,0.4,1,1,1];
in.mreinf   = ones(size(in.preinf,1),6);
in.spps     = (1/length(in.d))*ones(size(in.preinf,1),length(in.d));
in.ntrials  = 10000*ones(size(in.preinf,1),1);
in.ppun     = ones(size(in.preinf,1),6);
in.mpun     = in.ppun;

Val          = 3;    % Value of the integration window. The concept maps conceptually to the half decay time constant.
                     % Practical ex: If int. window is 5, then gamma =0.8, delta = 0.2.
                     % This concept provided me with a reasonable starting point to simulate
                     % algorithmic hypotheses that could go in line with the observed effects. 
Mod_factor   = 1/3;  % How much does the Laser affect either parameter. Value as in Lak et al., Neuron 2020
Gamma_delta  = 1;    % Correlation delta (1-gamma). Similar to fitted (see S3.4 in my thesis)

%Baseline(Control)
in.Intwindow = Val;
in.gamma     = (1-(1/ (in.Intwindow*Gamma_delta) ));          
in.delta     = (1/in.Intwindow);in.epsilon = 0.04;in.upsilon = (1/in.Intwindow);in.alpha = 0.1;in.zeta = 0.1;

for iSim = 1:nr_sim
[in,out(iSim)] = msimCDExperiment(in);
[DataMatrix] =Create_KDB_matrix(out(iSim));
nrStim = max(DataMatrix(:,1));
[Regressors,DependentVar, betas_B, pvalues, SEWeights, VIF] = fitGLM_Synth_datasets(DataMatrix,OneHotEncoder,TrialsInPast);
FittedStimulusMeans_B(iSim,:)       = betas_B(2:nrStim+1);
RewardedWeights_B(iSim,:)           = betas_B(nrStim+2:(nrStim+2+length(TrialsInPast)-1) );  %After rewarded trial
UnrewardedWeights_B(iSim,:)         = betas_B(nrStim+2+length(TrialsInPast):(nrStim+2+length(TrialsInPast)+length(TrialsInPast)-1)); %After unrewarded trial
end
std_B =std(RewardedWeights_B,1);
mean_B = mean(RewardedWeights_B,1);

%Hypothesis1: Reduced Integration Window (simultaneous increased delta and
%reduced gamma)
in.Intwindow   = Val-(Mod_factor)*Val;
in.gamma    = (1-(1/(in.Intwindow*Gamma_delta)));            
in.delta    = (1/in.Intwindow);in.epsilon = 0.04;in.upsilon = (1/in.Intwindow);in.alpha = 0.1;in.zeta = 0.1;

for iSim = 1:nr_sim
[in,out(iSim)] = msimCDExperiment(in);
[DataMatrix] =Create_KDB_matrix(out(iSim));
nrStim = max(DataMatrix(:,1));
[Regressors,DependentVar, betas_IW, pvalues, SEWeights, VIF] = fitGLM_Synth_datasets(DataMatrix,OneHotEncoder,TrialsInPast);
FittedStimulusMeans_IW(iSim,:)       = betas_IW(2:nrStim+1);
RewardedWeights_IW(iSim,:)           = betas_IW(nrStim+2:(nrStim+2+length(TrialsInPast)-1) );  %After rewarded trial
UnrewardedWeights_IW(iSim,:)         = betas_IW(nrStim+2+length(TrialsInPast):(nrStim+2+length(TrialsInPast)+length(TrialsInPast)-1)); %After unrewarded trial
end
std_IW =std(RewardedWeights_IW,1);
mean_IW = mean(RewardedWeights_IW,1);

%Hypothesis2: Increased delta
 in.gamma    = (1- 1/(Val*Gamma_delta));            
 in.delta    = 1/(Val-(Mod_factor)*Val);  %in.epsilon = 0.04;in.upsilon = 0.04;in.alpha = 0.1;in.zeta = 0.1;

for iSim = 1:nr_sim
[in,out(iSim)] = msimCDExperiment(in);
[DataMatrix] =Create_KDB_matrix(out(iSim));
nrStim = max(DataMatrix(:,1));
[Regressors,DependentVar, betas_D, pvalues, SEWeights, VIF] = fitGLM_Synth_datasets(DataMatrix,OneHotEncoder,TrialsInPast);
FittedStimulusMeans_D(iSim,:)       = betas_D(2:nrStim+1);
RewardedWeights_D(iSim,:)           = betas_D(nrStim+2:(nrStim+2+length(TrialsInPast)-1) );  %After rewarded trial
UnrewardedWeights_D(iSim,:)         = betas_D(nrStim+2+length(TrialsInPast):(nrStim+2+length(TrialsInPast)+length(TrialsInPast)-1)); %After unrewarded trial
end
std_D =std(RewardedWeights_D,1);
mean_D = mean(RewardedWeights_D,1);
%Hypothesis3: Increased gamma
 in.gamma    = (1- 1/( (Val*Gamma_delta) + (Mod_factor)* (Val*Gamma_delta) ));            
 in.delta    = 1/Val;

for iSim = 1:nr_sim
[in,out(iSim)] = msimCDExperiment(in);
[DataMatrix] =Create_KDB_matrix(out(iSim));
nrStim = max(DataMatrix(:,1));
[Regressors,DependentVar, betas_G, pvalues, SEWeights, VIF] = fitGLM_Synth_datasets(DataMatrix,OneHotEncoder,TrialsInPast);
FittedStimulusMeans_G(iSim,:)       = betas_G(2:nrStim+1);
RewardedWeights_G(iSim,:)           = betas_G(nrStim+2:(nrStim+2+length(TrialsInPast)-1) );  %After rewarded trial
UnrewardedWeights_G(iSim,:)         = betas_G(nrStim+2+length(TrialsInPast):(nrStim+2+length(TrialsInPast)+length(TrialsInPast)-1)); %After unrewarded trial
end
std_G =std(RewardedWeights_G,1);
mean_G = mean(RewardedWeights_G,1);
%% PLot

    figure;
    xaxis = TrialsInPast;
    ax1 = subplot(1,3,1);
    plot(xaxis,mean_B,'k-','LineWidth',2);hold on; plot(xaxis,mean_IW,'-','color',[0.9290 0.6940 0.1250],'LineWidth',2);
    
    %PLot shading STD baselines
    [x2,inBetween]= compute_shading_STD(mean_B, std_B);
    fill(x2, inBetween, 'k','FaceAlpha',0.3,'LineStyle','none');
    hold on;
    %PLot shading STD IW
    [x2,inBetween]= compute_shading_STD(mean_IW, std_IW);
    fill(x2, inBetween,[0.9290 0.6940 0.1250],'FaceAlpha',0.3,'LineStyle','none');
    
    axis(ax1,'square')
    axis([1-0.5 max(xaxis)+0.5 0.1 0.8 ])
    xlabel('Trials in the past')
    ylabel('Weights')
    title('Reduced Time Constant (Int. Win)')

    ax2 = subplot(1,3,2);
    plot(xaxis,mean_B,'k-','LineWidth',2);hold on; plot(xaxis,mean_D,'-','color',[0.9290 0.6940 0.1250]	,'LineWidth',2);
    %PLot shading STD baselines
    [x2,inBetween]= compute_shading_STD(mean_B, std_B);
    fill(x2, inBetween, 'k','FaceAlpha',0.3,'LineStyle','none');
    hold on;
    %PLot shading STD Delta
    [x2,inBetween]= compute_shading_STD(mean_D, std_D);
    fill(x2, inBetween,[0.9290 0.6940 0.1250],'FaceAlpha',0.3,'LineStyle','none');
    axis(ax2,'square')
    axis([1-0.5 max(xaxis)+0.5 0.1 0.8 ])
    title('Increase Delta')

    ax3 = subplot(1,3,3);
    plot(xaxis,mean_B,'k-','LineWidth',2);hold on; plot(xaxis,mean_G,'-','color',[0.9290 0.6940 0.1250]	,'LineWidth',2);
    %PLot shading STD baselines
    [x2,inBetween]= compute_shading_STD(mean_B, std_B);
    fill(x2, inBetween, 'k','FaceAlpha',0.3,'LineStyle','none');
    hold on;
    %PLot shading STD Gamma
    [x2,inBetween]= compute_shading_STD(mean_G, std_G);
    fill(x2, inBetween,[0.9290 0.6940 0.1250],'FaceAlpha',0.3,'LineStyle','none');
    axis(ax3,'square')
    axis([1-0.5 max(xaxis)+0.5 0.1 0.8 ])
    title('Increase Gamma')


%% HELPER FUNCTIONS
%% function [DataMatrix] =Create_KDB_matrix(out)
function [DataMatrix] =Create_KDB_matrix(out)

% Build DataMatrix
% 1 Stim seq
% 2 Cat 
% 3 Response (1 or 2)
% 4 Reward   (0 or 1)
% 5 Punishement (e.g, electric shock, reward omission does not count)

Correct_trials = sum(out.correct,2);
Potentially_reinf_trials = sum(out.reinf,2);
Reinf_trials   = Correct_trials.*Potentially_reinf_trials;
Session_vector = zeros(length(out.stimseq),1);
for iCond = 1:max(out.condition_vector)
    for jTrial = 1:length(out.stimseq)
        if out.condition_vector(jTrial) ==1 || out.condition_vector(jTrial)==2
            Session_vector(jTrial) = 1;
        elseif out.condition_vector(jTrial) ==3 || out.condition_vector(jTrial)==4
            Session_vector(jTrial) = 2;
        elseif out.condition_vector(jTrial) ==5 || out.condition_vector(jTrial)==6
            Session_vector(jTrial) = 3;
        end
    end
end

DataMatrix = zeros(length(out.stimseq),7);
DataMatrix(:,1) = out.stimseq;
DataMatrix(:,2) = out.category;
DataMatrix(:,3) = out.choices;
DataMatrix(:,4) = Reinf_trials;
DataMatrix(:,7) = out.condition_vector;
DataMatrix(:,6) = out.condition_vector;

end
%% function[x2,inBetween]= compute_shading_STD(mean, std_dev)
function[x2,inBetween]= compute_shading_STD(mean, std_dev)
    y = mean;
    x = 1:numel(mean);
    curve1 = y + std_dev;
    curve2 = y - std_dev;
    x2 = [x, fliplr(x)];
    inBetween = [curve1, fliplr(curve2)];
end
