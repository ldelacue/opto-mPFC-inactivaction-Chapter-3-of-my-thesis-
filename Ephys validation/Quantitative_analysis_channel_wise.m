%% Ephy analysis of optical inhibition of channel-wise mPFC activity

% Quantify effect of laser power (use 6.5 mW)
% Quantify effect of pulse duration (200,500,1000), cond 2, 3 & 4
% Quantify  extent of inhibition (channels).

% Needs mnspx.m, msep.m

%                                           Luis de la Cuesta 02/2025
%% Select recordings according to a selection criteria

load("All_curated_recordings.mat")
channels = zeros(5,50,length(myrecs));
Laser_power = 2; % 1 is 111, 2 is 206 and 3 is 400mW/mm2 upon which the selection is based
artOffset = 0.003; % in seconds, shift considered window by this size to avoid photoelectric artifacts.
for iRecs = 1:length(myrecs)
    myspx = myrecs{1,iRecs};
    mymat = nan(numel(myspx),2501,6);
    StimDur = [0.1,0.2,0.5,1,0.2,0.2]; % in s
    nrUnits(iRecs) = numel(myspx);
    for iUnits = 1:numel(myspx)
        UnitsInChannel(iUnits)= myspx(iUnits).channel;
        for jCond = 1:6
            mytrials_log = myspx(iUnits).Laserstims==myspx(iUnits).LaserInt_mW(Laser_power) & myspx(iUnits).codestims==jCond; % select appropriate Laser Power (6.5mW) & stim length
            trigtimes  = myspx(iUnits).timestims(mytrials_log)+artOffset;                % in s
            B_trigtimes  = myspx(iUnits).timestims(mytrials_log)+artOffset-StimDur(jCond); % in s

            spxtimes   = myspx(iUnits).timespx;
            pre = 0;
            post = (StimDur(jCond)-artOffset)*1000;% in ms due to requirements of mnspx
            nspx = mnspx(spxtimes,trigtimes,pre,post);
            b_nspx = mnspx(spxtimes,B_trigtimes,pre,post);

            [h(iUnits,jCond),p(iUnits,jCond)] = ttest (nspx, b_nspx);
            FR_Base(iUnits,jCond) = sum(b_nspx)/((StimDur(jCond)-2*artOffset)*(length(b_nspx)));   %in s
            FR_Stim(iUnits,jCond) = sum(nspx)/((StimDur(jCond)-2*artOffset)*(length(nspx)));       %in s
        end
    end
    FR_Base_Nor = FR_Base./FR_Base;
    FR_Stim_Nor = FR_Stim./FR_Base;
    Modulation = sum(h(:,1:4),2);
    Modulated_id =(Modulation>=3); %chose those units that at least have been statistically modulated in three (out of 4) conditions.
    % figure; plot([ones(length(FR_Base_Nor(Modulated_id)),1),FR_Stim_Nor(Modulated_id,1:4)]')
    UnitsInChannel(Modulated_id)
    UnitsInChannel(~Modulated_id);
    channels(1,1:length(unique(UnitsInChannel)),iRecs) = unique(UnitsInChannel);
    for iChan = 1:length(channels)
        channels(2,iChan,iRecs) = length(find(UnitsInChannel==channels(1,iChan,iRecs)))    ;            % number of units per channel
        channels(3,iChan,iRecs) = length(find(UnitsInChannel(Modulated_id)==channels(1,iChan,iRecs)));  % number of modulated units per channel
    end
    channels(4,:,iRecs) = channels(3,:,iRecs)./channels(2,:,iRecs);                                       % percentage of modulated units per channel
    channels(5,:,iRecs) = channels(4,:,iRecs)>=0.5;                                                       % Modulated(1) or not (0)

    if sum(channels(5,:,iRecs)) >=4 % if the recording has less than 4 significant modulated channels do not consider the recording
        Keep_recordings(1,iRecs) = sum(channels(5,:,iRecs)); %display score
        Keep_recordings(2,iRecs) = 1;                        %display outcome
        disp('Keep this recording')
    else
        disp('Discard this recording')
        Keep_recordings(1,iRecs) = sum(channels(5,:,iRecs)); % display score
        Keep_recordings(2,iRecs) = 0;                        % display outcome
    end
    clearvars -except myrecs channels Keep_recordings iRecs artOffset Laser_power nrUnits
end

%% categorical scatter plot
Keeps = channels(:,:,Keep_recordings(2,:)==1);
iKeeps = 1;
iChann = 2;
channel_modulation = nan(32,size(Keeps,3));
for iKeeps= 1:size(Keeps,3)
    for iChann=1:32
        try
        channel_modulation(iChann,iKeeps) = Keeps(5,Keeps(1,:,iKeeps)==iChann,iKeeps);
        catch
        end
    end
end
channel_modulation = flip(channel_modulation);
channel_modulation =channel_modulation(1:20,:)*1;

% Plot
figure
for iChann = 1:size(channel_modulation,2)
Chann1 = channel_modulation(channel_modulation(:,1)==1,iChann);
x = iChann*ones(1, length(channel_modulation(:,iChann)));
y = (1:1:length(channel_modulation(:,iChann)))*-100;
Mod_chann_log = channel_modulation(:,iChann)==1;
Not_Mod_chann_log = channel_modulation(:,iChann)==0;
scatter(x(Mod_chann_log)', y(Mod_chann_log)',100,[0.9290 0.6940 0.1250],'filled');hold on
scatter(x(Not_Mod_chann_log), y(Not_Mod_chann_log),100, [0.4940 0.1840 0.5560], 'filled'); hold on
end
hold on ; line([0 5],[0 0],'Color','black','LineStyle','--')
axis ([0.5, 4.5, (length(channel_modulation)*-100)-100, 0])
xticks([1 2 3 4])
xticklabels({'LF64 Lh P','LF64 Rh P','LF65 Rh A','LF66 Lh'})
ylabel('Distance to optic fiber (um)')
axis square
title('Modulated channels (500ms & 206 mW/mm^2)')
%% Analyize the data to generate recordings-wise inactivation plots
Valid_recordings = find((Keep_recordings(2,:)==1));
MeanFRs = [];
SemFRs  = [];
Modulation_rates = [];
SEP_Mod_rates =    [];
for ivRecs = 1:length(Valid_recordings)
    myspx = myrecs{1,Valid_recordings(ivRecs)};
    StimDur = [0.1,0.2,0.5,1,0.2,0.2]; % in s
    for iLInt  = 1:3
        for iUnits = 1:numel(myspx)
            KeepUnits(iUnits) = myspx(iUnits).channel>12;
            UnitsInChannel(iUnits)= myspx(iUnits).channel;
            for iCond = 1:6
                mytrials_log = myspx(iUnits).Laserstims==myspx(iUnits).LaserInt_mW(iLInt) & myspx(iUnits).codestims==iCond; % select appropriate Laser Power (6.5mW) & stim length
                trigtimes  = myspx(iUnits).timestims(mytrials_log)+artOffset;                % in s
                B_trigtimes  = myspx(iUnits).timestims(mytrials_log)+artOffset-StimDur(iCond); % in s
                spxtimes   = myspx(iUnits).timespx;
                pre = 0;
                post = (StimDur(iCond)-artOffset)*1000;% in ms due to requirements of mnspx
                nspx = mnspx(spxtimes,trigtimes,pre,post);
                b_nspx = mnspx(spxtimes,B_trigtimes,pre,post);
                [h(iUnits,iCond),p(iUnits,iCond)] = ttest (nspx, b_nspx);
                FR_Base(iUnits,iCond) = sum(b_nspx)/((StimDur(iCond)-2*artOffset)*(length(b_nspx)));   %in s
                FR_Stim(iUnits,iCond) = sum(nspx)/((StimDur(iCond)-2*artOffset)*(length(nspx)));       %in s
            end
        end
        FR_Stim_Nor = FR_Stim./FR_Base;

        % Analysis 1: How strongly are modulated units modulated per stim length and laser power?
        Bad_apples = sum( isnan(h(:,1:4)),2);
        Bad_apples = Bad_apples~=0;
        Modulation = sum(h(:,1:4),2);
        Modulation = Modulation(~Bad_apples,:);
        h          = h(~Bad_apples,1:4);
        Modulated_id =h==1; %chose those units that at least have been statistically modulated in three (out of 4) conditions.
        FR_Stim_Nor = FR_Stim_Nor(~Bad_apples,1:4);
        KeepUnits    = KeepUnits(~Bad_apples);


        % restrict which units you use to the channels that are at least
        % from 2mm from the optic fiber (selected in KeepUnits).
        FR_Stim_Nor = FR_Stim_Nor(KeepUnits,:);
        Modulated_id= Modulated_id(KeepUnits,:);
        Modulation  = Modulation(KeepUnits);
        h           = h(KeepUnits,:);

        for iCond = 1:4
            if isnan(mean( FR_Stim_Nor( Modulated_id(:,iCond),iCond),1))
                MeanFR(1,iCond)  = 0;
            else
                MeanFR(1,iCond) =mean( FR_Stim_Nor( Modulated_id(:,iCond),iCond),1);
            end
            SemFR(1,iCond)   = std ( FR_Stim_Nor( Modulated_id(:,iCond),iCond),1) / sqrt( size (FR_Stim_Nor( Modulated_id(:,iCond),iCond),1) ) ;
        end
        MeanFRs  = vertcat(MeanFRs, MeanFR ); %vertically  concatenate results from different inhibition powers (4,4)
        SemFRs   = vertcat(SemFRs,  SemFR );

        % Analysis 2: How many modulated units are there per stim length and laser power?

        Modulation_rate   = sum(h(:,1:4),1)/length(Modulation);
        SEP_Mod_rate      = msep(Modulation_rate,length(Modulation)*ones(1,length(Modulation_rate)));
        Modulation_rates  = vertcat(Modulation_rates,Modulation_rate);
        SEP_Mod_rates     = vertcat(SEP_Mod_rates,SEP_Mod_rate);
    end
    clearvars -except Valid_recordings MeanFRs SemFRs Modulation_rates SEP_Mod_rates ivRecs myrecs artOffset
end

%% PLotting
xaxis= [1:1:5];
figure
subplot(2,2,1)
for ivRecs = 1:length(Valid_recordings)
    errorbar(xaxis,[1,MeanFRs(4+ivRecs,:)],[0,SemFRs(4+ivRecs,:)]) % this selects all rows with the power of 206mW (rows 5:9)
    valuesForRegression(ivRecs,:) = MeanFRs(4+ivRecs,:);
    hold on
end
xticks([1 2 3 4 5])
axis square
xticklabels({'B','100','200','500','1000'})
xlabel('Inhibition window (ms)')
ylabel('Norm. firing rate')
title('Laser Power = 206 mW/mm^2')
axis([0.5 max(xaxis)+0.5 0 1])

subplot(2,2,2)
MeanFRs_pow = [MeanFRs(1:4,3),MeanFRs(5:8,3),MeanFRs(9:12,3)]; % we select the 3rd window length(500ms), because is the closest to the actual parameter used for behavioural testing
SemFRs_pow  = [SemFRs(1:4,3),SemFRs(5:8,3),SemFRs(9:12,3)];
xaxis= [1:1:4];
for ivRecs = 1:length(Valid_recordings)
    errorbar(xaxis,[1,MeanFRs_pow(ivRecs,:)],[0,SemFRs_pow(ivRecs,:)])
    valuesForRegression_pow(ivRecs,:) = MeanFRs_pow(ivRecs,:);

    hold on
end
xticks([1 2 3 4])
axis square
xticklabels({'B','111','206','400'})
xlabel('Power (mw/mm2)')
ylabel('Norm. firing rate')
title('Inhibition window = 500 ms')
axis([0.5 max(xaxis)+0.5 0 1])

subplot(2,2,3)
xaxis= [1:1:5];
for ivRecs = 1:length(Valid_recordings)
    errorbar(xaxis,[0,Modulation_rates(4+ivRecs,:)],[0,SEP_Mod_rates(4+ivRecs,:)])   
    valuesForRegression_MR(ivRecs,:) = Modulation_rates(4+ivRecs,:);
    hold on
end
xticks([1 2 3 4 5])
axis square
xticklabels({'B','100','200','500','1000'})
xlabel('Inhibition window (ms)')
ylabel('fraction of modulated units')
title('Laser Power = 206 mW/mm^2')
axis([0.5 max(xaxis)+0.5 0 1])

subplot(2,2,4)
Modulation_rates_pow = [Modulation_rates(1:4,3),Modulation_rates(5:8,3),Modulation_rates(9:12,3)]; % we select the 3rd window length(500ms), because is the closest to the actual parameter used for behavioural testing
SEP_Mod_rates_pow  = [SEP_Mod_rates(1:4,3),SEP_Mod_rates(5:8,3),SEP_Mod_rates(9:12,3)];
xaxis= [1:1:4];
for ivRecs = 1:length(Valid_recordings)
    errorbar(xaxis,[0,Modulation_rates_pow(ivRecs,:)],[0,SEP_Mod_rates_pow(ivRecs,:)])
    valuesForRegression_MR_pow(ivRecs,:) = Modulation_rates_pow(ivRecs,:);
    hold on
end
xticks([1 2 3 4])
axis square
xticklabels({'B','111','206','400'})
xlabel('Power (mw/mm2)')
ylabel('fraction of modulated units')
title('Inhibition window = 500 ms')
axis([0.5 max(xaxis)+0.5 0 1])
legend({'LF64 Lh P','LF64 Rh P','LF65 Rh A','LF66 Lh'})
