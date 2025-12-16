function [in,out] = msimCDExperiment(uin)
% function [in,out] = msimCDExperiment(uin)
%
% Simulates a behavioral experiment for the KDB model, the Lak model, the Lak model with an additional punishment parameter,
% a randomly drifting criterion, and a fixed criterion.
% For the drifting criterion, we add a random number drawn from a normal distribution with µ=0 and sd = in.delta on each trial, but the criterion is
% decremented by gamma (if not, it will run off).
% The fixed criterion model corresponds to classic signal detection theory.
%
% Allows to specify different contingencies with varying stimulus presentation probabilities (SPPs), reward probabilities and magnitudes, and
% punishment probabilities and magnitudes.
%
% Plot the results with plotCDSim(out).
%
% INPUTS
% All inputs parameters should be fields of the struct 'uin'
% d             vector with means of stimulus distributions, default [-1,1]
% category      vector with category assignment of stimuli, default [1,2] - must have the same number of elements as d
% c             starting value of criterion, relevant only for models 'drift' and 'fixed', default 0
% model         which model to use to update the criterion; default 'income'
%               'income', 'unrewarded', 'error', 'general' (income + unrewarded), 'Lak', 'LakPun', 'drift', 'fixed'
% ntrials       number of trials per condition, default 10000
%               if preinf and ppun have several rows, ntrials must have the same number of rows,
%               where each row specifies the number of trials per condition
% spps          vector specifying stimulus presentation probabilities, default: uniform (must add to 1 of course)
% preinf        vector specifying reinforcer probabilities for each stimulus, default [1,1]; as in in.spps, each row denotes a different condition
% mreinf        vector specifying reinforcer magnitude, default [1,1]
% mpun          vector specifying punisher magnitude, default [1,1]
% ppun          vector specifying punishment probabilities for each stimulus, default [1,1]; each row denotes a different condition
%               note: punishment is given for incorrect responses only and mpun and ppun need to be >0 for the error-based model to work
%
% gamma         leak rate, default 0.99
% delta         learning rate for reinforcers, default 0.04
% epsilon       learning rate for errors or punishers, default 0.04
% upsilon       learning rate for unrewarded trials, default 0.04
% alpha         learning rate parameter for the Lak model, default 0.1
% zeta          learning rate parameter for the LakPun model, default 0.1
% decnoise      parameter for decision noise, default 0
%
% OUTPUT
% in            struct with all parameters (defaults and from uin)
% out           struct with fields 'stimseq', 'c' (criterion values),'choices','category','correct','reinf' (reinforcement or not per trial and category),
%               'pun' (see 'reinf'), 'confidence', 'condition_vector', 'blocks' (for plotting)
%
% Maik C. Stüttgen, March 2024, built on simkdbsession3.m
%% specify default settings and overwrite if necessary
in.d = [-1,1];in.category = [1,2];in.c = 0;in.model ='income';in.nconditions = 1;in.ntrials = 1000;
in.spps = ones(1,numel(in.d))/numel(in.d);in.preinf = [1,1];in.mreinf = [1,1];
in.ppun = [1,1];in.mpun = [1,1];
% Int Window
% in.Intwindow   = 3.5;
% in.gamma    = (1-(1/in.Intwindow)); % 0.98;
% in.delta    = (1/in.Intwindow);in.epsilon = 0.04;in.upsilon = (1/in.Intwindow);in.alpha = 0.1;in.zeta = 0.1;


% Delta
%  in.gamma    = 1- 1/3.5;
%  in.delta    = 1/2.5;  %in.epsilon = 0.04;in.upsilon = 0.04;in.alpha = 0.1;in.zeta = 0.1;

% Gamma
%  in.gamma    = 1- 1/4.5;
%  in.delta    = 1/3.5;

in.decnoise = 0;
if ~exist('uin','var'),uin = in;disp('No input provided, running code with default values.'),end
in = overwriteDefaults(in,uin); % also does some input checks
clear uin
if strcmp(in.model,'income SLR RD')
    in.delta = in.delta.*[0.1,0.3,1,1,0.3,0.1];
    in.delta = in.delta*1.5;
elseif  strcmp(in.model,'income SLR')
    in.delta = in.delta.*[0.1,0.3,1,1,0.3,0.1];
    in.delta = in.delta*1.5;
end
%% inform user about what you are going to do
disp('------------------------------------------------------------------')
if strcmp(in.model,'income'),disp(['SIMULATING INCOME-BASED MODEL WITH ',num2str(sum(in.ntrials)),' TRIALS'])
elseif strcmp(in.model,'unrewarded'),disp(['SIMULATING UNREWARDED-BASED MODEL WITH ',num2str(sum(in.ntrials)),' TRIALS'])
elseif strcmp(in.model,'error'),disp(['SIMULATING ERROR-BASED MODEL WITH ',num2str(sum(in.ntrials)),' TRIALS'])
elseif strcmp(in.model,'general'),disp(['SIMULATING GENERAL MODEL WITH ',num2str(sum(in.ntrials)),' TRIALS'])
elseif strcmp(in.model,'Lak'),disp(['SIMULATING THE LAK MODEL WITH ',num2str(sum(in.ntrials)),' TRIALS'])
elseif strcmp(in.model,'LakPun'),disp(['SIMULATING THE LAK MODEL WITH PUNISHMENT WITH ',num2str(sum(in.ntrials)),' TRIALS'])
elseif strcmp(in.model,'drift'),disp(['SIMULATING A RANDOMLY DRIFTING CRITERION MODEL WITH ',num2str(sum(in.ntrials)),' TRIALS'])
elseif strcmp(in.model,'fixed'),disp(['SIMULATING A MODEL WITH A FIXED CRITERION WITH ',num2str(sum(in.ntrials)),' TRIALS'])
end

%% generate stimulus sequence
breaks        = cumsum(in.ntrials);
thisCondition = 1;
for i = 1:in.nconditions

    % adjust (increase) number of trials in this condition so that all n stimuli can be presented according to their pre-specified probabilities
    % first, check whether an integer number of trials is produced in this condition
    % if not, increase number of trials so that it is
    if ~all(rem(in.spps(i,:)*in.ntrials(i),1)==0)
        in.ntrials(i) = ceil(in.ntrials(i)/(numel(in.d)*40))*(numel(in.d)*40);
        disp(['Increasing in.ntrials(i) to ',num2str(sum(ceil(in.spps(2,:)*in.ntrials(i))))])

        % adjust variable breaks to new number of trials in this condition
        if i==1,breaks(i) = in.ntrials(i);
        elseif i>1,breaks(i) = in.ntrials(i)+breaks(i-1);end

    end

    % generate stimulus sequence
    allStims = [];
    for j = 1:numel(in.d)
        allStims = [allStims;j*ones(in.ntrials(i)*in.spps(i,j),1)];
    end
    out.stimseq((i-1)*in.ntrials(i)+1:breaks(i),1) = shuffle(allStims);
end
clear allStims

%% prepare and preallocate output field arrays
out.c          = zeros(sum(in.ntrials)+1,1);
out.choices    = nan(sum(in.ntrials),1);
out.category   = nan(sum(in.ntrials),1);      % needed for kdbfit
out.correct    = zeros(sum(in.ntrials),1);    % needed for kdbfit
out.reinf      = zeros(sum(in.ntrials),2);    % for Rf1 and Rf2
out.pun        = zeros(sum(in.ntrials),2);    % for Pun1 and Pun2
out.confidence = zeros(sum(in.ntrials),1);

% for speed, generate vectors of evidence values, category, and decision noise here already
out.EVs      = normrnd(in.d(out.stimseq),1)';
out.category = in.category(out.stimseq)';
out.decnoise = in.decnoise*randn(sum(in.ntrials),1);

% needed for the Lak model only
out.V_R1 = 0.5*ones(sum(in.ntrials),1);   % action value R1; V_0 is set to 0.5 for optimistic reward expectation
out.V_R2 = 0.5*ones(sum(in.ntrials),1);   % action value R2; V_0 is set to 0.5 for optimistic reward expectation
out.RPE  = nan(sum(in.ntrials),1);        % reward prediction error (delta)
out.negRPE = nan(sum(in.ntrials),1);      % negative reward prediction error (for LakPun)
out.Q_C  = nan(sum(in.ntrials),1);        % chosen values
out.p_C1 = 1-normcdf(out.EVs,0,1);        % confidence category 1 for all trials
out.p_C2 = 1-out.p_C1;                    % confidence category 2 for all trials

%% run the simulation
for i = 1:sum(in.ntrials)

    % which condition are we in?
    if i>breaks(thisCondition)
        thisCondition = thisCondition + 1;
    end

    % first, determine value of DV on current trial, add decision noise and compute decision confidence
    currentDV = out.EVs(i);
    out.c(i) = out.c(i)+out.decnoise(i);
    out.confidence(i,1) = currentDV - out.c(i);

    % second, determine choice of subject on this trial
    if ~strcmp(in.model,'Lak') && ~strcmp(in.model,'LakPun') && currentDV<out.c(i)
        out.choices(i) = 1;

    elseif ~strcmp(in.model,'Lak') && ~strcmp(in.model,'LakPun') && currentDV>out.c(i)
        out.choices(i) = 2;

    elseif strcmp(in.model,'Lak') || strcmp(in.model,'LakPun')
        Q_R1 = out.V_R1(i) * out.p_C1(i);
        Q_R2 = out.V_R2(i) * out.p_C2(i);
        if Q_R2>Q_R1
            out.choices(i,1) = 2;
            out.Q_C(i,1) = Q_R2;
        else
            out.choices(i,1) = 1;
            out.Q_C(i,1) = Q_R1;
        end

    end

    % third, note down consequences (out.reinf, out.pun, and out.correct are initialized with 0s)
    if in.category(out.stimseq(i))==out.choices(i) % correct R1 or R2
        out.correct(i,out.choices(i)) = 1;
        if rand<in.preinf(thisCondition,out.stimseq(i))
            if out.choices(i)==1
                out.reinf(i,1) = in.mreinf(thisCondition,1);
            elseif out.choices(i)==2
                out.reinf(i,2) = in.mreinf(thisCondition,2);
            end
        end
    elseif in.category(out.stimseq(i))~=out.choices(i) % incorrect R1 or R2
        if rand<in.ppun(thisCondition,out.stimseq(i))
            if out.choices(i)==1
                out.pun(i,1) = in.mpun(thisCondition,1);
            elseif out.choices(i)==2
                out.pun(i,2) = in.mpun(thisCondition,2);
            end
        end
    end

    % last, update decision criterion
    if strcmp(in.model,'income')  % original
        out.c(i+1) = in.gamma*out.c(i) + in.delta*(out.reinf(i,1)-out.reinf(i,2));

    elseif strcmp(in.model,'income RD')  % adjusted by Luis to include no update no pullback
        if out.reinf(i,1)==1 || out.reinf(i,2)==1
            out.c(i+1) = in.gamma*out.c(i) + in.delta*(out.reinf(i,1)-out.reinf(i,2));
        else
            out.c(i+1) = out.c(i);
        end
    elseif strcmp(in.model,'income SLR')  % original
        out.c(i+1) = in.gamma*out.c(i) + in.delta(out.stimseq(i))*(out.reinf(i,1)-out.reinf(i,2));
    elseif strcmp(in.model,'income SLR RD')  % adjusted by Luis to include no update no pullback
        if out.reinf(i,1)==1 || out.reinf(i,2)==1
            out.c(i+1) = in.gamma*out.c(i) + in.delta(out.stimseq(i))*(out.reinf(i,1)-out.reinf(i,2));
        else
            out.c(i+1) = out.c(i);
        end

    elseif strcmp(in.model,'unrewarded')
        if any(out.reinf(i,:)==1)
            out.c(i+1) = in.gamma*out.c(i);
        else
            if out.choices(i,1)==1
                out.c(i+1) = in.gamma*out.c(i) + in.upsilon*-1;
            elseif out.choices(i,1)==2
                out.c(i+1) = in.gamma*out.c(i) + in.upsilon;
            end
        end

    elseif strcmp(in.model,'error')
        out.c(i+1) = in.gamma*out.c(i) + in.epsilon*(out.pun(i,2)-out.pun(i,1));

    elseif strcmp(in.model,'general')
        out.c(i+1) = in.gamma*out.c(i) + in.delta*(out.reinf(i,1)-out.reinf(i,2)) ...
            + in.upsilon*(out.pun(i,2)-out.pun(i,1));

    elseif strcmp(in.model,'Lak')
        out.RPE(i) = sum(out.reinf(i,:)) - out.Q_C(i);
        if out.choices(i,1)==1
            out.V_R1(i+1) = out.V_R1(i) + in.alpha*out.RPE(i);
            out.V_R2(i+1) = out.V_R2(i);
        elseif out.choices(i,1)==2
            out.V_R1(i+1) = out.V_R1(i);
            out.V_R2(i+1) = out.V_R2(i) + in.alpha*out.RPE(i);
        end

    elseif strcmp(in.model,'LakPun')
        out.RPE(i) = sum(out.reinf(i,:)) - out.Q_C(i);
        out.negRPE(i) = out.Q_C(i) - sum(out.pun(i,:));
        if out.choices(i)==1
            out.V_R1(i+1) = out.V_R1(i) + in.alpha*out.RPE(i) - in.zeta*out.negRPE(i);
            out.V_R2(i+1) = out.V_R2(i);
        elseif out.choices(i)==2
            out.V_R1(i+1) = out.V_R1(i);
            out.V_R2(i+1) = out.V_R2(i) + in.alpha*out.RPE(i) - in.zeta*out.negRPE(i);
        end

    elseif strcmp(in.model,'drift')
        out.c(i+1) = in.gamma*out.c(i) + randn*in.delta;

    elseif strcmp(in.model,'fixed')
        out.c(i+1) = out.c(i);

    end

end

out.condition_vector = ones(size(out.stimseq,1),1);
cumulated_trials = cumsum(in.ntrials);
for i = 1:in.nconditions
    out.condition_vector((i-1)*in.ntrials(i)+1:cumulated_trials(i)) = i;
end
clear cumulated_trials i

end
%% NESTED FUNCTION OVERWRITEDEFAULTS
function in = overwriteDefaults(in,uin)
uinfieldnames = fieldnames(uin);
for i = 1:numel(uinfieldnames)
    switch uinfieldnames{i}
        case 'd',in.d = uin.d;
        case 'category',in.category = uin.category;
        case 'c',in.c = uin.c;
        case 'model',in.model = uin.model;
        case 'ntrials',in.ntrials = uin.ntrials;
        case 'spps',in.spps = uin.spps;
        case 'preinf',in.preinf = uin.preinf;
        case 'mreinf',in.mreinf = uin.mreinf;
        case 'ppun',in.ppun = uin.ppun;
        case 'mpun',in.mpun = uin.mpun;
        case 'gamma',in.gamma = uin.gamma;
        case 'delta',in.delta = uin.delta;
        case 'upsilon',in.upsilon = uin.upsilon;
        case 'epsilon',in.epsilon = uin.epsilon;
        case 'alpha',in.alpha = uin.alpha;
        case 'zeta',in.zeta = uin.zeta;
        case 'decnoise',in.decnoise = uin.decnoise;
    end
end
in.nconditions = size(in.spps,1);

infieldnames = fieldnames(in);
ndims = nan(numel(infieldnames),2);
for i = 1:numel(infieldnames)
    ndims(i,:) = size(in.(infieldnames{i}));  % a nice way to avoid eval
end

if numel(unique(ndims(6:11,1)))>1
    disp([num2str(ndims(:,1),'%4.0f_'),num2str(ndims(:,2),'%4.0f_ '),char(infieldnames)])
    error('The input fields do not specify the same number of conditions.')
elseif numel(unique(ndims([1,2,7,8,9,10,11],2)))>1
    disp([num2str(ndims(:,1),'%4.0f_'),num2str(ndims(:,2),'%4.0f_ '),char(infieldnames)])
    error('The input fields do not specify the same number of stimuli.')
elseif any(sum(in.spps,2)<0.9999) || any(sum(in.spps,2)>1.0001)
    disp(num2str(sum(in.spps,2)))
    error('Stimulus presentation probabilities do not sum to 1.')
end

end
