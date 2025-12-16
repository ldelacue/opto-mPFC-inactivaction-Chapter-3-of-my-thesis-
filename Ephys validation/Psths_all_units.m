%% PSTHs ALL UNITS IN A RECORDING
% PLot all psths of all recording units in each recording in each of the 6 different conditions
% (100ms laser, 200ms, 500ms, 1000ms, LED, LED + 200ms laser).
% Based on a maik analysis code.

% Needs mpsth

                                                % Luis de la Cuesta 02/2025
clear all
close all
clc

load('All_curated_recordings.mat')
myspx = LF64_Lh_P_struct; %choose one of the recordings, Lh stands for left hemis. and P for posterior recording.
mymat = nan(numel(myspx),2501,6);

for iUnits = 1:numel(myspx)
    for jCond = 1:6
        psth = mpsth(myspx(iUnits).timespx,myspx(iUnits).timestims(myspx(iUnits).codestims==jCond),'fr',1,'tb',0,'pre',1000,'post',1500);
        mymat(iUnits,:,jCond) = msdf(psth,'exp',10,'notext');
    end
end
%%
% mean
figure,hold on
titles = { '100ms', '200ms', '500ms', '1000ms', 'LED', 'LED + 200ms'};

for iCond = 1:6
    subplot(2,3,iCond),hold on
    plot(-1000:1500,mean(mymat(:,:,iCond)))
    plot([0,0],[0,8],'k:')
    axis([-1000,1000,0,8])
    title(titles{iCond})
end

%% example unit
for iUnit = 1:size(mymat,1)
Nr_unit = iUnit;
titles = { '100ms', '200ms', '500ms', '1000ms', 'LED', 'LED + 200ms'};
Opto_Shading_x = [100,200,500,1000,NaN,200]; % End x value

for iCond = 1:6
    
    figure('name',['peri-stimulus time histogram (Unit Nr ',num2str(iUnit),')'],'units','normalized','position',[0.3 0.4 0.4 0.2])
    [psth,trialspx] = mpsth_opto(myspx(Nr_unit).timespx,myspx(Nr_unit).timestims(iCond:6:end),'pre',500,'post',1200,'binsz',20,'chart',2);
    title(titles{iCond})
  % Opto_Shading_y = [min(psth(:,1)), max(psth(:,1))+10];
    Raster_Shading_y = [length(myspx(Nr_unit).timestims(iCond:6:end))+10, length(myspx(Nr_unit).timestims(iCond:6:end))+10];
    subplot(211)
    area( [0,Opto_Shading_x(iCond)]', Raster_Shading_y','FaceColor','y','FaceAlpha',.2)
    subplot(212)
    hold on
     Histogram_Shading_y = [max(psth(:,2))+5, max(psth(:,2))+5];
     area( [0,Opto_Shading_x(iCond)]', Histogram_Shading_y','FaceColor','y','FaceAlpha',.2)

end 

end


%% HELPER FUNCTIONS

%% function [psth,trialspx] = mpsth_opto(spxtimes,trigtimes,varargin)
function [psth,trialspx] = mpsth_opto(spxtimes,trigtimes,varargin)

% [psth,trialspx] = mpsth(spxtimes,trigtimes,varargin)
% function generates a peri-stimulus time histogram (psth) with time base in column 1 and histogram in column 2
% in addition, function returns spike timestamps relative to trigger times
% IMPORTANT: all timestamp inputs (spxtimes, trigtimes) must be seconds and will be converted to ms in the script!
% IMPORTANT: the psth will extend from -pre:binsz:post+binsz!!!
% 
% MANDATORY INPUTS
% spxtimes      vector with timestamps (seconds) of spike events
% trigtimes     vector with timestamps (seconds) of trigger events
%
% OPTIONAL INPUTS
% pre           time before trigger to include in psth (default 1000 ms)
% post          time after trigger to include in psth (default 1000 ms)
% fr            if '1', normalizes to firing rate (Hz); if '0', does nothing (default)
% tb            if '1', function returns timebase in the first column (default); if '0', no time base - output is a single column
% binsz         bin size of psth (default: 1 ms)
%               note that the frequencies/spike counts in a psth will include all spikes within binsz ms after the respective time point
%               thus, when binsz is 10 and psth(x,1)=-500, the value in psth(x,2) corresponds to the time range -500 to -491 ms pre-event
% chart         if '0' (default), no plot will be generated
%               if '1', a PSTH will be generated
%               if '2', a PSTH together with a raster plot will be generated
%
% EXAMPLES
% load unitForMLIBTesting.mat
% psth = mpsth(spx.timings,spx.eventtimings(spx.eventmarkers==9),'fr',1,'chart',2);
% 
% psth = mpsth(chan9.timings(chan9.markers(:,1)==1),chan32.timings(chan32.markers(:,1)==105))
%               generates a psth (time base in first column, psth in second column) from marker 1 in channel 9 for event 105
%
% [psth trialspx] = mpsth(chan9.timings(chan9.markers(:,1)==1),chan32.timings(chan32.markers(:,1)==105),'pre',3000,'post',3000,'fr',1);
%               same, but PSTH extends from -3 to +3 s around event (rather than -1 to +1 s, which is the default)
%               scales to firing rate ('fr' set to 1)
%
% HISTORY
% feb 1, 2019    tried out another algorithm for psth contruction, compared speed to the other two
%                went back to do the for-loop because the code from two days ago was slower
%                also, I debugged the program for binsz>1, there were issues with the very last time bin
% jan 29, 2019   moved 'remove time base' to the end of the code because charts could not be plotted when tb was 0
%                faster code for generating psths (single line rather than for-loop)
% jan 28, 2019   new code for generating psth; legacy code still contained
%                added input check for binsz to "define and override defaults", removed it from "construct psth & trialspx new..."
%                replaced try with "if ~isempty(trialspx{i})" in "construct psth & trialspx new..."
% dec 4, 2015    ylabel for PSTH now matches user request for firing rate or spike count
% mar 31, 2015   mraster eliminated from code because of inefficiency and imprecision
% feb 25, 2015   minor bug detected: when trigtimes contains nans, firing rates are not scaled appropriately -> fixed
% feb 19, 2015   trialspx now holds spike time stamps with precision > binsz; until today, spikes times in trialspx were rounded and thus
%                allocated to the nearest integer value (with binsz=1, 0.9 would end up as 1, 0.2 as 0)
%                now each bin contains spikes from t to t+bin rather than t-bin/2 to t+bin/2
% feb 12, 2015   removed +binsz from psth plot
% sep 12, 2013   minor bug detected; in preallocation, it previously read "psth = zeros(pre/binsz+post/binsz,2)", which
%                could result in mismatches of matrix size
% july 30, 2013  minor changes to comments above
% feb 23, 2012   changed conversion of spike counts into firing rate -> bug eliminated
% august 12      changed argument name 'ntb' to 'tb' for consistency
% april 18       debugged bin size argument, added optional er display
% april 16       added bin size as argument
% feb 16, 2011   added examples
%
% by Maik C. Stuettgen, Feb 2011
%% define and override defaults
spxtimes  = spxtimes*1000;
trigtimes = trigtimes(~isnan(trigtimes))*1000;
pre   = 1000;
post  = 1000;
fr    = 1;
tb    = 1;
binsz = 1;
chart = 0;
Opto_Shading_x = [100,200,500,1000,NaN,200]; % End x value

if nargin
  for i=1:2:size(varargin,2)
    switch varargin{i}
      case 'pre'
        pre = varargin{i+1};
      case 'post'
        post = varargin{i+1};
      case 'fr'
        fr = varargin{i+1};
      case 'tb'
        tb = varargin{i+1};
      case 'binsz'
        binsz = varargin{i+1};
      case 'chart'
        chart = varargin{i+1};
      otherwise
        errordlg('unknown argument')
    end
  end
end
if binsz<1 || binsz>1000 || mod(binsz,1)~=0
  error('Invalid bin size specified');
end
%% pre-allocate for speed
if binsz>1
  if rem(pre,binsz)~=0 || rem(post,binsz)~=0
    error('Both pre-stim and post-stim windows for the PSTH must be divisible by binsz without remainder.')
    % Otherwise, none of the windows would extend from 0 to binsz.
  else
    psth(:,1) = (-1*pre:binsz:post); % time base
    psth(:,2) = 0;
  end
elseif binsz==1
  % in this case, pre+post+1 bins are generated (ranging from pre:1:post)
  psth = zeros(pre+post+1,2);
  psth (:,1) = (-1*pre:1:post);       % time base
end
trialspx = cell(numel(trigtimes),1);
%% construct psth & trialspx new - used from January 28, 2019 - is much faster than old code
% rationale: when we have -5:5 as time, we have 11 bins in which to count spikes
% we count spikes >=pre and <post+1 so we get the exact number of bins
% if we had <=post+1, we would have to have one more bin when the timestamp=post
% for construction of the psth, every timestamp in trialspx is then assigned its 'floored' value
% psth2 = psth;
% psth3 = psth;
for i = 1:numel(trigtimes)
  trialspx{i} = spxtimes((spxtimes - trigtimes(i))>=-pre & (spxtimes - trigtimes(i))<(post+binsz))-trigtimes(i);
  if ~isempty(trialspx{i})
    % this is the old way to do things, and it turned out to be fastest
    for j = 1:numel(trialspx{i})
      psth(floor(trialspx{i}(j)/binsz+pre/binsz+1),2) = psth(floor(trialspx{i}(j)/binsz+pre/binsz+1),2)+1;
    end
    % alternative code, slowest
%     % now need to check whether any indices appear twice; if so, we need to run the for-loop...
%     if numel(floor(trialspx{i}/binsz+pre/binsz+1))==numel(unique(floor(trialspx{i}/binsz+pre/binsz+1)))
%       psth(floor(trialspx{i}/binsz+pre/binsz+1),2) = psth(floor(trialspx{i}/binsz+pre/binsz+1),2)+1;  % this line is new, the for-loop is legacy!
%     else
%       for j = 1:numel(trialspx{i})
%         psth(floor(trialspx{i}(j)/binsz+pre/binsz+1),2) = psth(floor(trialspx{i}(j)/binsz+pre/binsz+1),2)+1;
%       end
%     end
  end
end
% allspx = cell2vectors(trialspx);
% psth3(:,2) = hist(allspx(:,1),-pre+binsz/2:binsz:post+binsz/2);     % this is intermediate
% figure,hold on
% plot(psth(:,1),psth(:,2)),plot(psth2(:,1),psth2(:,2)),plot(psth3(:,1),psth3(:,2))
%% construct psth & trialspx legacy - USED UNTIL JANUARY 28, 2019
% for i = 1:numel(trigtimes)
%   spikes = spxtimes - trigtimes(i);                    % all spikes relative to current trigtime
%   trialspx{i} = spikes(spikes>=-pre & spikes<=post+1);   % spikes close to current trigtime
%   if binsz>0  % this is strange - binsz must be 1 or larger
%     try       % why try?
%       for j = 1:numel(trialspx{i})
%         psth(floor(trialspx{i}(j)/binsz+pre/binsz+1),2) = psth(floor(trialspx{i}(j)/binsz+pre/binsz+1),2)+1;
%       end
%     end
%   else
%     error('Invalid bin size specified');
%   end
%   clear spikes
% end
%% normalize to firing rate if desired
if fr==1
  psth (:,2) = (1/binsz)*1000*psth(:,2)/numel(trigtimes);
end
%% plot
if chart==1
  %figure('name',['peri-stimulus time histogram (bin size ',num2str(binsz),' ms)'],'units','normalized','position',[0.3 0.4 0.4 0.2])
  bar(psth(:,1),psth(:,2),'k','BarWidth',1)
  axis([min(psth(:,1))-10 max(psth(:,1))+10 0 max(psth(:,2))+1])
  xlabel('peri-stimulus time (ms)')
  if fr
    ylabel('spikes per s')
  else
    ylabel('counts per bin')
  end
elseif chart==2
  %figure('name',['peri-stimulus time histogram (bin size
  %',num2str(binsz),' ms)'],'units','normalized','position',[0.3 0.3 0.4 %0.3])                                   % Luis edit to plot a different title
  subplot(212)
  bar(psth(:,1)+binsz/2,psth(:,2),'k','BarWidth',1)
  axis([min(psth(:,1))-10 max(psth(:,1))+10 0 max(psth(:,2))+1])
  %axis([min(psth(:,1))-10 max(psth(:,1))+10 0 70+1])

  
  xlabel('peri-stimulus time (ms)')
  if fr==1
    ylabel('spikes per s')
  else
    ylabel('counts per bin')
  end
  subplot(211),hold on  % make rastter
  for i = 1:numel(trialspx)
    plot(trialspx{i},i*ones(numel(trialspx{i}),1),'k.','MarkerSize',2,'LineStyle','none')
  end
  axis([-pre-10 post+10 0 numel(trialspx)+1])
  ylabel('trials')


end
%% remove time base
if tb==0
  psth(:,1) = [];
end
end
