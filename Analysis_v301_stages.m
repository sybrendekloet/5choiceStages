function Analysis_v301_stages
%% Main script for analysis of 5-CSRTT learning stage data
% 
% Ideally, this script will be all you need to analyze the entire dataset
% of a 5-CSRTT learning stages experiment.
%
% It will only require you to load the appropriate data file, which is
% generally called something like: "20200609_allData.mat"
%
% You can then run every chapter to get behavioral and neurophysiological
% results.

%% Table of contents
% 1) Pre-processing
%   1.1) Initiation and housekeeping
%   1.2) Data organization
% 2) Behavior
% 3) Photometry


%% 1.1) Initiation and housekeeping
% Define key variables
placements = {'GDM' 'GVM'}; % indicators for different experimental groups
conditionNames = {'stage1' 'stage2' 'stage3' 'stage4' 'stage5' ...
    'stage6' 'stage7' 'stage8' 'fixedITI' 'varSD' 'varITI'}; % indicators for different stages of task

% Where figures of analysis will be saved
exptname = 'LearningStages'; % name folder to which experiment analysis will belong
savefolder = ['D:\Data\FPCompiled_stages\Analysis\' exptname]; % folder where this analysis will be saved
mkdir(savefolder); % creates directory where figures will end up



%% 1.2) Data organization
% Generate and organize traces from pre-processed data file

% Will take the loaded datafile as input (input variable is named
% 'compiledData'). Loops through every rat that has data in the file, and
% for every session, the script will cut the session calcium trace into
% smaller traces around task-relevant moments (i.e. trials).
% 
% Output:
% tmpTrc    - matrix with signal traces of each trial from all rats
% ana_meta  - matrix with decoding parameters for each trial


% Define key variables
rNms = fieldnames(compiledData); 
ratNames = rNms(~contains(rNms,'placements')); % Generate cell of all rat names in experiment
field_order = {'GDM' 'GVM'}; % Define order of groups (required for later steps)
tmpTrc = cell(1,2); baseTrc = cell(1,2); corrTmpTrc = cell(1,2); ana_meta = zeros(1,11); % Initiate a bunch of arrays that will be filled in this chapter
sync_win = {ceil(-5*15.89):ceil(12.5*15.89), ...
    ceil(-12.5*15.89):ceil(5*15.89)}; % Define trace length. Holds two windows, one for trials synced at trial start, and one for trials synced at response
trialN = 40; % Threshold number of trials
accuracyT = .8; % Threshold accuracy
omissionT = .2; % Threshold omissions
errC = 0; % counter for sessions with few baseline traces
spec = [5 7.5 12.5;...
    1 0.5 0.2]; % Indicators for ITI times (1st row) and SD times (2nd row)
tic % Start timer

for rat = 1:numel(ratNames) % Loops through all rats in dataset
    
    toc % Report start time of loop for each rat. If this takes excessively long, contact Syb.    
    group = find(contains(field_order, placements(strcmp(ratNames{rat}(isletter(ratNames{rat})),placements)))==1); %Determine group rat belongs to    
    conditions = fieldnames(compiledData.(ratNames{rat})); %Define of which stages we have data in this rat
    
    for cond = 1:numel(conditions) % Loops through all stages this rat has done
        
        sessN = numel(compiledData.(ratNames{rat}).(conditions{cond}).msfin); % Define number of sessions rat has done in each stage        
        condCheck = find(strcmp(conditions{cond},conditionNames)==1); % Match stage with appropriate stage number
        trained = zeros(1,1); % Reset 'trained' variable to 0
        
        for sess = 1:sessN % Loops through all sessions rat has done in current stage
            
            if isnan(compiledData.(ratNames{rat}).(conditions{cond}).trialstart{sess}) % Check if session actually has trials, if not skip
                continue
            elseif size(compiledData.(ratNames{rat}).(conditions{cond}).trialstart{sess},1)<5 % Also skip if less than 5 trials in this session (data not reliable)
                continue                
            else  
                % Generate traces
                data = compiledData.(ratNames{rat}).(conditions{cond}).msfin{sess}; % Select data
                meta = [compiledData.(ratNames{rat}).(conditions{cond}).trialstart{sess}(:,1),...
                    compiledData.(ratNames{rat}).(conditions{cond}).response{sess}(:,1),...
                    compiledData.(ratNames{rat}).(conditions{cond}).trialstart{sess}(:,2:4)]; % Generate meta-matrix that holds all trial parameters like response type, or latency
                
                % Exclude trials too close to session start or end
                validtr = find(meta(:,1)>12.5 & meta(:,2)<(numel(data)/15.89)-20);
                
                % Index parameter for data matrix size. This variable holds info on size of the
                % matrix, so we can easily append data to the bottom after every iteration
                % of the loop.
                msize = size(tmpTrc{1},1);
                
                % Sometimes it can be helpful to normalize the data. These
                % variables store the 2nd and 98th percentile of dF/F trace
                % of the current session, to store it into the meta-matrix
                % for later use if desired.
                upper_prc = prctile(data,98);
                lower_prc = prctile(data,2);
                
                % Assigns variable ITI and variable SD trials to an ITI
                % duration or cue length.
                % 1: easy trials (5s ITI or 1s SD; baseline)
                % 2: medium trials (7.5s ITI or 0.5s SD)
                % 3: hard trials (12.5s ITI or 0.2s SD)
                if mean(meta(validtr,4))>0
                    if mean(meta(validtr,4))>1 % Tests if session is vITI 
                        [~,subspec]=ismember(meta(validtr,4),spec(1,:));
                    else % Tests if session is vSD
                        [~,subspec]=ismember(meta(validtr,4),spec(2,:));
                    end
                else % Adds a zero if 
                    subspec = zeros(numel(validtr),1);
                end
                
                % Trialstart and response times relative to session length
                trtime = meta(validtr,1);
                rsptime = meta(validtr,2);
                
                for sync = 1:2 % Loop through different synchronization moments (1 = trialstart, 2 = response)
                    % Generate trial traces
                    iidx = bsxfun(@plus, ceil(meta(validtr,sync)*(1017.3/64)), sync_win{sync}); % Create all peri-trial indices
                    tmpTrc{sync}(msize+1:msize+numel(validtr),:) = data(iidx); % Generate traces
                    
                    % Naive/non-trial related signal
                    if condCheck<2 % Check if stage 1, this requires an exception because trialstart and response are the same
                        iidx2 = bsxfun(@plus, ceil(meta(validtr,sync)*(1017.3/64)), ceil(-2*15.89):ceil(2*15.89)); % Create all indices
                        c = [[1;iidx2(:,end)],[iidx2(:,1);numel(data)]]; % Non-trial window edges. ...
                        % This is aimed to also store data points that are not near any trials, to look at naive baseline behavior.
                    else
                        c = [[1;iidx(:,end)],[iidx(:,1);numel(data)]]; %non-trial window edges
                    end
                    
                    % Create a matrix with baseline calcium signals. These
                    % signals are not related to any trial dependent
                    % activity.
                    swin = ceil(5*15.89); % Window size
                    ttmp = zeros(size(iidx,1),swin); % Initiate temporary variable that stores baseline traces
                    for id = 1:size(c,1)-1 % Loops through non-trial windows during session
                        d = c(id,1):c(id,2); % Difference between end of previous trial and start of next
                        if (c(id,2)-c(id,1))>swin % Test if time between two trials is enough to take a baseline window
                            ws = randi([d(1), d(end)-ceil(swin)],1); % Randomly pick start of baseline window from period between two trials
                            ttmp(id,:) = data(ws:ws-1+swin); % Generate trace
                        else
                            ttmp(id,:) = nan(1,swin); % Add NaN to matrix if window is too short for baseline trace
                        end
                    end
                    baseTrc{sync}(msize+1:msize+numel(validtr),:)=ttmp; % Store baseline calcium signal in matrix
                    if sum(~isnan(ttmp(:,1)))<5 % If there are fewer than 5 baseline traces, report error and continue
                        errC = errC+1
                    end
                    corrTmpTrc{sync}(msize+1:msize+numel(validtr),:)=(data(iidx)-nanmean(nanmean(ttmp,2)))./nanstd(nanmean(ttmp,2)); % Store baseline-corrected trial traces
                    % Final Trace = (calcium traces)-mean(sessionbaseline)
                    
%                     % Plots
%                     if (sync == 1 && cond <10)
%                         figure(rat)
%                         set(gcf,'position',get(0,'screensize'));
%                         subplot(3,3,cond)
%                         plot(movmean(mean(data(iidx)),25),'color',col_rep(sess),'linewidth',2)
%                         hold on
%                         meanE = mean(mean(ttmp(~isnan(ttmp(:,1)),:)));
%                         stE = std(mean(ttmp(~isnan(ttmp(:,1)),:)));
%                         line([1,size(data(iidx),2)],[meanE meanE], 'linestyle', '--', 'color', col_rep(sess));
%                         jbfill(1:size(data(iidx),2),ones(1,size(data(iidx),2))*(meanE-stE), ...
%                             ones(1,size(data(iidx),2))*(meanE+stE), col_rep(sess), col_rep(sess),1,0.2);
%                         hold on
%                         if sess == sessN
%                             yl = get(gca, 'ylim');
%                             line([5*15.89 5*15.89],yl, 'color','k')
%                             line([10*15.89 10*15.89],yl, 'color','k')
%                             set(gca,'ylim',yl)
%                             xlim([1 12.5*15.89])
%                         end
%                     end
                    
                end % sync
                
                % Primer for distinction between trained or untrained
                % trials performed in each stage. See below for
                % explanation ("Explanation A")
                if cond < 4 
                    if numel(validtr)<trialN
                        trained = zeros(numel(validtr),1);
                    else
                        trained = [zeros(trialN,1); ones(numel(validtr)-trialN,1)];
                    end
                else 
                    trained = ones(numel(validtr),1)*-1;
                end
                
                % Here, the meta-matrix is constructed. Columns correspond
                % to:
                % 1      2    3     4     5     6        7           8
                % group, rat, cond, sess, resp, subspec, lower norm, upper norm
                % 9             % 10        11
                % trial time,   resp time   trained/untrained
                ana_meta(msize+1:msize+numel(validtr),:) = ...
                    [ones(numel(validtr),1)*group,...
                    ones(numel(validtr),1)*rat, ones(numel(validtr),1)*condCheck,...
                    ones(numel(validtr),1)*sess, meta(validtr,3),...
                    subspec,...
                    ones(numel(validtr),1)*lower_prc, ones(numel(validtr),1)*upper_prc,...
                    trtime,rsptime,trained];
            end    
        end
        
        % Explanation A: 
        % For each stage above stage 2, assign training status (trained/untrained)
        % and append to trial data matrix (ana_meta)
        % Definition of 'trained':
        %      > 80% accuracy, <20 omissions, 40 trials per stage
        % The script will move a sliding window (40 trials) across the stage and look
        % for the trial where performance criteria are met. All trials
        % before that point will be labeled 'untrained' (0 in column 11 in meta matrix). Trials beyond that
        % point are labeled 'trained' (1 in column 11 in meta matrix).
        
        if condCheck > 3 % Check if stage is above stage 3. (Before that, we can't really apply the performance thresholds)
            trid_in = (ana_meta(:,1)==group & ana_meta(:,2)==rat & ana_meta(:,3)==condCheck); % ID of trials in current stage
            trials_in = ana_meta(trid_in,:); % Trial metadata of trials in current stage
            
            if condCheck == 9 && sum(ana_meta(:,1)==group & ana_meta(:,2)==rat & ana_meta(:,3)<condCheck)==0
                % Check if rat has done all training stages before (some have no training data). If not -> stage 9 only has 'trained' trials
                trained = ones(sum(ana_meta(:,1)==group & ana_meta(:,2)==rat & ana_meta(:,3)==condCheck),1);
                    ana_meta(trid_in,end)=trained;
            
            elseif size(trials_in,1)>5 % Check if more than 5 trials in current stage (otherwise continue; not reliable)
                accuracy = 0; omissions = 1; tr =1; % Set start values for sliding window.
                
                while (accuracy<accuracyT || omissions > omissionT) && tr<(size(trials_in,1)-trialN)
                    accuracy = sum(trials_in(tr:tr+trialN,5)==1)/sum(trials_in(tr:tr+trialN,5)<3); % Sliding window for accuracy %
                    omissions = sum(trials_in(tr:tr+trialN,5)==3)/trialN; % Sliding window for omission %
                    tr = tr+1; % Trial counter
                end
                
                if (tr-1+trialN)==size(trials_in,1) % If performance threshold was never met, put '0'
                    ana_meta(trid_in,end)=zeros(size(trials_in,1),1);
                else
                    
                    % Append meta matrix with training status
                    trained = [zeros(tr-1+trialN,1); ones(size(trials_in,1)-(tr-1+trialN),1)];
                    ana_meta(trid_in,end)=trained;
                end
            end
        end
    end % Stage/cond
end % Rat

data_in{1} = tmpTrc{1}; % Store data in different variable (to make sure it isn't overwritten when running again)
data_in{2} = tmpTrc{2};
meta_in = ana_meta; % Store metadata in different variable
[meta_in_sorted, id_order] = sortrows(sortrows(meta_in,3),2);
data_in_sorted{1} = data_in{1}(id_order);
data_in_sorted{2} = data_in{2}(id_order);

% Z-Transform the data (on single trial level)
bwin = [ceil(1*15.89)^0, ceil(4*15.89)]; % baseline window = -5s : -1s from trialstart
data_in_z{1} = (data_in{1}-mean(data_in{1}(:,bwin(1):bwin(2)),2))./...
    std(data_in{1}(:,bwin(1):bwin(2)),[],2); % Z-correction for trialstart-synced traces
data_in_z{2} = (data_in{2}-mean(data_in{1}(:,bwin(1):bwin(2)),2))./...
    std(data_in{1}(:,bwin(1):bwin(2)),[],2); % Z-correction for response-synced traces


%% 2) Behavior
% This chapter will analyze behavioral data to answer basic questions about
% 5-CSRTT learning performance. 
%
% Input: meta matrix called 'meta_in', generated in previous chapter
%
% Output: 
%   2.0) Structure array called 'behavior'. Holds response parameters, 
%        latencies, post-error effects
%        Each field holds 3 matrices: the first two correspond to dmPFC rats and
%        vmPFC rats, respectively. The third holds data from all rats.
% 
% 6 figures of different behavioral analyses
%   Part 1: Just learning stages
%   2.1) How many trials does it take to reach stage performance
%   thresholds?
%   2.2) Raw numbers of each response type per stage
%   2.3) Indirect behavior parameters (1)
%   2.4) Indirect behavior parameters (2)
%   2.5a)Sliding window of behavioral performance across stages
%   2.5b)Sliding window across stages, trained and untrained trials
%   separated
%   2.5c)Difference between trained/untrained within and across stages
%
%   Part 2: variable ITI and SD included
%
%
%   2.6a) Predictive power of learning stages in sessions that require higher
%        cognitive load (1)
%   2.6b) Predictive power of learning stages in sessions that require higher
%        cognitive load (2)
%
% Other useful things:
% Meta table contents (useful when tweaking this chapter)
% Generate or append trial data matrix
% 1      2    3     4     5     6        7           8
% group, rat, cond, sess, resp, subspec, lower norm, upper norm
% 9             % 10         % 11
% trial time,   resp time    learning/nolearning

behav_in = meta_in_sorted;
ncond = max(unique(behav_in(:,3))); % Number of distinct stages in dataset
behavior = struct; % Create structure array for behavioral data
delays = [5 7.5 12.5]; % vITI delay times

for gr = 1:3
    if gr<3
        nsize = unique(behav_in(behav_in(:,1)==gr,2)); % id of all rats in group
    else
        nsize = unique(behav_in(behav_in(:,1)<=2,2)); % id of all rats in group
    end
    for ra = 1:numel(nsize)
        for st = 1:ncond
            if gr<3
            trid2 =  (behav_in(:,1)==gr & behav_in(:,2)==nsize(ra) & ...
                behav_in(:,3)==st);                
            else
            trid2 =  (behav_in(:,1)<=2 & behav_in(:,2)==nsize(ra) & ...
                behav_in(:,3)==st);
            end
            if sum(trid2)==0
                continue
            end
            
            meta2 = behav_in(trid2,:);
            k = 1;
            % Store values
            % Means
            if st>2 && st<10
                %Raw numbers
                behavior.correct{gr}(st,ra,k)=sum(meta2(:,5)==1);
                behavior.incorrect{gr}(st,ra,k)=sum(meta2(:,5)==2);
                behavior.omission{gr}(st,ra,k)=sum(meta2(:,5)==3);
                behavior.premature{gr}(st,ra,k)=sum(meta2(:,5)==4);
                
                %With defined performance threshold
                behavior.trialsBeforeTrained{gr}(ra,st)=sum(behav_in(trid2,end)==0);
                behavior.trialsAfterTrained{gr}(ra,st)=sum(behav_in(trid2,end)==1);
                
                % Proportional numbers
                behavior.accuracy{gr}(st,ra,k)=sum(meta2(:,5)==1)/(sum(meta2(:,5)<3));
                behavior.omissionpct{gr}(st,ra,k)=sum(meta2(:,5)==3)/(size(meta2,1));
                behavior.prematurepct{gr}(st,ra,k)=sum(meta2(:,5)==4)/(size(meta2,1));
               
                % Latencies
                behavior.correctlatency{gr}(st,ra,k)=mean(meta2(meta2(:,5)==1,10)-meta2(meta2(:,5)==1,9)-5);
                behavior.incorrectlatency{gr}(st,ra,k)=mean(meta2(meta2(:,5)==2,10)-meta2(meta2(:,5)==2,9)-5);
                behavior.prematuretime{gr}(st,ra,k)=mean(meta2(meta2(:,5)==4,10)-meta2(meta2(:,5)==4,9));
                
                % Post error latencies
                meta_err = [0;meta2(1:end-1,5)];
                behavior.posterrorlatency{gr}(st,ra,k)=mean(meta2(meta_err>1 & meta2(:,5)==1,10)-meta2(meta_err>1 & meta2(:,5)==1,9)-5);
                
            elseif st<=2
                behavior.accuracy{gr}(st,ra,k)=sum(meta2(:,5)==1)/(sum(meta2(:,5)<3));
                behavior.correct{gr}(st,ra,k)=sum(meta2(:,5)==1);
                behavior.correctlatency{gr}(st,ra,k)=mean(meta2(meta2(:,5)==1,10)-meta2(meta2(:,5)==1,9));
                behavior.trialsBeforeTrained{gr}(ra,st)=sum(meta_in(trid2,end)==0);
                behavior.trialsAfterTrained{gr}(ra,st)=sum(meta_in(trid2,end)==1);
            
%             elseif st==1
%                 behavior.accuracy{gr}(st,ra,k)=sum(meta2(:,5)==1)/(sum(meta2(:,5)<3));
%                 behavior.correct{gr}(st,ra,k)=sum(meta2(:,5)==1);
%                 behavior.trialsBeforeTrained{gr}(ra,st)=sum(meta_in(trid2,end)==0);
%                 behavior.trialsAfterTrained{gr}(ra,st)=sum(meta_in(trid2,end)==1);
                
            elseif st>9
                for k = 1:3
                    % Store vSD and vITI parameters (in 3D matrix for
                    % different categories)
                    behavior.correct{gr}(st,ra,k)=sum(meta2(:,5)==1 & meta2(:,6)==k);
                    behavior.incorrect{gr}(st,ra,k)=sum(meta2(:,5)==2 & meta2(:,6)==k);
                    behavior.omission{gr}(st,ra,k)=sum(meta2(:,5)==3 & meta2(:,6)==k);
                    behavior.premature{gr}(st,ra,k)=sum(meta2(:,5)==4 & meta2(:,6)==k);
                    
                    behavior.accuracy{gr}(st,ra,k)=sum(meta2(:,5)==1 & meta2(:,6)==k)/(sum(meta2(:,5)<3 & meta2(:,6)==k));
                    behavior.omissionpct{gr}(st,ra,k)=sum(meta2(:,5)==3 & meta2(:,6)==k)/(size(meta2(meta2(:,6)==k),1));
                    behavior.prematurepct{gr}(st,ra,k)=sum(meta2(:,5)==4 & meta2(:,6)==k)/(size(meta2(meta2(:,6)==k),1));
                
                    if st==11 %vITI
                        behavior.correctlatency{gr}(st,ra,k)=mean(meta2(meta2(:,5)==1 & meta2(:,6)==k,10)-meta2(meta2(:,5)==1 & meta2(:,6)==k,9)-delays(k));
                        behavior.incorrectlatency{gr}(st,ra,k)=mean(meta2(meta2(:,5)==2 & meta2(:,6)==k,10)-meta2(meta2(:,5)==2 & meta2(:,6)==k,9)-delays(k));
                        behavior.prematuretime{gr}(st,ra,k)=mean(meta2(meta2(:,5)==4 & meta2(:,6)==k,10)-meta2(meta2(:,5)==4 & meta2(:,6)==k,9));
                        % Post error
                        meta_err = [0;meta2(1:end-1,5)];
                        behavior.posterrorlatency{gr}(st,ra,k)=mean(meta2(meta_err>1 & meta2(:,5)==1 & meta2(:,6)==k,10)-meta2(meta_err>1 & meta2(:,5)==1 & meta2(:,6)==k,9)-delays(k));
                        
                    else %vSD
                        behavior.correctlatency{gr}(st,ra,k)=mean(meta2(meta2(:,5)==1 & meta2(:,6)==k,10)-meta2(meta2(:,5)==1 & meta2(:,6)==k,9)-5);
                        behavior.incorrectlatency{gr}(st,ra,k)=mean(meta2(meta2(:,5)==2 & meta2(:,6)==k,10)-meta2(meta2(:,5)==2 & meta2(:,6)==k,9)-5);
                        behavior.prematuretime{gr}(st,ra,k)=mean(meta2(meta2(:,5)==4 & meta2(:,6)==k,10)-meta2(meta2(:,5)==4 & meta2(:,6)==k,9));
                        
                        % Post error
                        meta_err = [0;meta2(1:end-1,5)];
                        behavior.posterrorlatency{gr}(st,ra,k)=mean(meta2(meta_err>1 & meta2(:,5)==1 & meta2(:,6)==k,10)-meta2(meta_err>1 & meta2(:,5)==1 & meta2(:,6)==k,9)-5);
                        
                    end
                end
            end            
        end
        
        % Store single trials
        if gr<3
            trid3{1} =  (meta_in_sorted(:,1)==gr & meta_in_sorted(:,2)==nsize(ra)); % All trials
            trid3{2} =  (meta_in_sorted(:,1)==gr & meta_in_sorted(:,2)==nsize(ra) & meta_in_sorted(:,11)==0); % Untrained trials
            trid3{3} =  (meta_in_sorted(:,1)==gr & meta_in_sorted(:,2)==nsize(ra) & meta_in_sorted(:,11)==1); % Trained trials
        else
            trid3{1} =  (meta_in_sorted(:,1)<=gr & meta_in_sorted(:,2)==nsize(ra)); % All trials
            trid3{2} =  (meta_in_sorted(:,1)<=gr & meta_in_sorted(:,2)==nsize(ra) & meta_in_sorted(:,11)==0); % Untrained trials
            trid3{3} =  (meta_in_sorted(:,1)<=gr & meta_in_sorted(:,2)==nsize(ra) & meta_in_sorted(:,11)==1); % Trained trials
        end
        behavior.singleTrials{gr}{ra,1}=meta_in_sorted(trid3{1},:); % Store all trials
        behavior.singleTrials{gr}{ra,2}=meta_in_sorted(trid3{2},:); % Store untrained trials
        behavior.singleTrials{gr}{ra,3}=meta_in_sorted(trid3{3},:); % Store trained trials
        
        % Sliding performance window for all trials
        ws = 0; maxws = 20;
        meta3 = meta_in_sorted(trid3{1},:);
        for tx = 1:sum(trid3{1})
            behavior.accwin{gr}{ra,1}(tx) = sum(meta3(tx-ws:tx,5)==1)/sum(meta3(tx-ws:tx,5)<=2);
            behavior.corwin{gr}{ra,1}(tx) =  sum(meta3(tx-ws:tx,5)==1)/sum(meta3(tx-ws:tx,5)<=4);
            behavior.incwin{gr}{ra,1}(tx) =  sum(meta3(tx-ws:tx,5)==2)/sum(meta3(tx-ws:tx,5)<=4);
            behavior.omwin{gr}{ra,1}(tx) =  sum(meta3(tx-ws:tx,5)==3)/sum(meta3(tx-ws:tx,5)<=4);
            behavior.prwin{gr}{ra,1}(tx) =  sum(meta3(tx-ws:tx,5)==4)/sum(meta3(tx-ws:tx,5)<=4);
            if ws<maxws
                ws =ws+1; 
            end
        end
        
        for trn = 1:2 % Store trials based on trained or untrained status
            behavior.accwin{gr}{ra,trn+1} = behavior.accwin{gr}{ra,1}(meta3(:,11)==(trn-1));
            behavior.corwin{gr}{ra,trn+1} =  behavior.corwin{gr}{ra,1}(meta3(:,11)==(trn-1));
            behavior.incwin{gr}{ra,trn+1} = behavior.incwin{gr}{ra,1}(meta3(:,11)==(trn-1));
            behavior.omwin{gr}{ra,trn+1} =  behavior.omwin{gr}{ra,1}(meta3(:,11)==(trn-1));
            behavior.prwin{gr}{ra,trn+1} = behavior.prwin{gr}{ra,1}(meta3(:,11)==(trn-1));            
        end
        
    end
    
    % Remove zeros
    behavior.trialsBeforeTrained{gr}(behavior.trialsBeforeTrained{gr}==0)=nan;
    behavior.trialsAfterTrained{gr}(behavior.trialsAfterTrained{gr}==0)=nan;
end

%% Figures (Only learning stages)
% Set folder to save figures
savefolder2 = 'D:\Publications\Learning paper\Raw Figures'; % folder where this analysis will be saved
mkdir(savefolder2); % creates directory where figures will end up
close all % Close all previously opened figures

% 1. How many trials does it take to reach threshold performance?
% a) During stages
figure(1)
hold on
% boxplot(behavior.trialsBeforeTrained{1},'positions',0.75:1:8.75, 'plotstyle','compact','colors','b','jitter',0.3);
% boxplot(behavior.trialsBeforeTrained{2},'positions',1:9, 'labels', conditionNames(1:9),'plotstyle','compact','colors','r','jitter',0.3);
boxplot(behavior.trialsBeforeTrained{3},'positions',1.25:1:9.25, 'labels', conditionNames(1:9),'plotstyle','compact','colors','k','jitter',0.3);
yl = get(gca, 'ylim');
line([8.5 8.5],yl, 'color','k','linestyle','-.')
ylim([0 yl(2)]);
labels = findobj(gca,'Tag','Box');
legend([labels(end), labels(ceil(numel(labels)/2)), labels(1)],'dmPFC', 'vmPFC', 'all', 'location', 'northwest')
ylabel('Trials')
set(gca, 'box','off')

% subplot(2,1,2)
% % With individual rats (Doesn't rly work due to missing data)
% plot(cumsum(behavior.trialsBeforeTrained{1},2)')


% 2. Raw numbers of each response type per stage
figure(2)
vnames = {'correct' 'incorrect' 'omission' 'premature'};
for pm = 1:numel(vnames)
    subplot(2,4,pm+4*(gr-1))
    behavior.(vnames{pm}){gr}(behavior.(vnames{pm}){gr}==0)=nan;
    hold on
    boxplot(behavior.(vnames{pm}){1}(1:9,:,1)','positions',0.75:1:8.75, 'labels', conditionNames(1:9),'plotstyle','compact','colors',col_rep(1));
    boxplot(behavior.(vnames{pm}){2}(1:9,:,1)','positions',1:9, 'labels', conditionNames(1:9),'plotstyle','compact','colors',col_rep(2));
    boxplot(behavior.(vnames{pm}){3}(1:9,:,1)','positions',1.25:1:9.25, 'labels', conditionNames(1:9),'plotstyle','compact','colors','k');
    yl(1,:) = get(gca,'ylim');
    
    %vSD
    atmp=behavior.(vnames{pm}){gr}(10,:,:);
    btmp=[atmp(:,:,1); atmp(:,:,2); atmp(:,:,3)];
    boxplot(btmp', 'positions', 9.75:.25:10.25,'plotstyle','compact','colors',col_resp(pm))
    yl(2,:) = get(gca,'ylim');
    
    %vITI
    ctmp=behavior.(vnames{pm}){gr}(11,:,:);
    dtmp=[ctmp(:,:,1); ctmp(:,:,2); ctmp(:,:,3)];
    boxplot(dtmp', 'positions', 10.75:.25:11.25,'plotstyle','compact','colors',col_resp(pm))
    xlim([0 12])
    
    yl(3,:) = get(gca, 'ylim');
    line([9.5 9.5],yl, 'color','k','linestyle','-.')
    line([10.5 10.5],yl, 'color','k','linestyle','-.')
    ylim([min(yl(:,1)), max(yl(:,2))]);
    ylabel('Trials')
    title(vnames{pm})
end

% 3. Indirect behavior parameters (1)
figure(3)
vnames = {'accuracy' 'omissionpct' 'prematurepct'};
for pm = 1:numel(vnames)
    subplot(numel(vnames),1,pm)
%     behavior.(vnames{pm}){1}(behavior.(vnames{pm}){1}==0)=nan;
%     behavior.(vnames{pm}){2}(behavior.(vnames{pm}){2}==0)=nan;
    behavior.(vnames{pm}){3}(behavior.(vnames{pm}){3}==0)=nan;
    hold on
%     boxplot(behavior.(vnames{pm}){1}(1:9,:,1)','positions',0.8:1:8.8,'labels', conditionNames(1:9), 'plotstyle','compact','colors',col_rep(1));
%     boxplot(behavior.(vnames{pm}){2}(1:9,:,1)','positions',1:9,'labels', conditionNames(1:9), 'plotstyle','compact','colors',col_rep(2));
    boxplot(behavior.(vnames{pm}){3}(1:9,:,1)','positions',1.2:1:9.2,'labels', conditionNames(1:9), 'plotstyle','compact','colors','k');
    yl23 = get(gca,'ylim');
    
    %vSD
    %         atmp=behavior.(vnames{pm}){gr}(10,:,:);
    %         btmp=[atmp(:,:,1); atmp(:,:,2); atmp(:,:,3)];
    %         boxplot(btmp', 'positions', 9.75:.25:10.25,'plotstyle','compact','colors',col_resp(pm))
    %         yl(2,:) = get(gca,'ylim');
    %
    %         %vITI
    %         ctmp=behavior.(vnames{pm}){gr}(11,:,:);
    %         dtmp=[ctmp(:,:,1); ctmp(:,:,2); ctmp(:,:,3)];
    %         boxplot(dtmp', 'positions', 10.75:.25:11.25,'plotstyle','compact','colors',col_resp(pm))
    %         xlim([0 12])
    
%     yl(3,:) = get(gca, 'ylim');
%     line([9.5 9.5],[0 1], 'color','k','linestyle','-.')
%     line([10.5 10.5],[0 1], 'color','k','linestyle','-.')
    ylim([min(yl23), 1.05*max(yl23)]);
    ylabel('Proportion')
    title(vnames{pm})
    set(gca,'box','off')

    if pm ==3
        labels = findobj(gca,'Tag','Box');
        legend([labels(end), labels(ceil(numel(labels)/2)), labels(1)],'dmPFC', 'vmPFC', 'all', 'location', 'northwest')
    end
end


% 4. Indirect behavior parameters (2)
figure(4)
vnames = {'correctlatency' 'incorrectlatency' 'posterrorlatency' 'prematuretime'};

for gr= 1:2
    for pm = 1:numel(vnames)
        subplot(2,4,pm+4*(gr-1))
        behavior.(vnames{pm}){gr}(behavior.(vnames{pm}){gr}==0)=nan;
        hold on
        boxplot(behavior.(vnames{pm}){1}(3:9,:,1)','positions',3:9,'labels', conditionNames(3:9), 'plotstyle','compact','colors',col_resp(pm));
        yl(1,:) = get(gca,'ylim');
        
         %vSD
        atmp=behavior.(vnames{pm}){gr}(10,:,:);
        btmp=[atmp(:,:,1); atmp(:,:,2); atmp(:,:,3)];
        boxplot(btmp', 'positions', 9.75:.25:10.25,'plotstyle','compact','colors',col_resp(pm))
        yl(2,:) = get(gca,'ylim');
        
        %vITI
        ctmp=behavior.(vnames{pm}){gr}(11,:,:);
        dtmp=[ctmp(:,:,1); ctmp(:,:,2); ctmp(:,:,3)];
        boxplot(dtmp', 'positions', 10.75:.25:11.25,'plotstyle','compact','colors',col_resp(pm))
        xlim([0 12])
        
        yl(3,:) = get(gca, 'ylim');
        line([9.5 9.5],yl, 'color','k','linestyle','-.')
        line([10.5 10.5],yl, 'color','k','linestyle','-.')
        ylim([min(yl(:,1)), max(yl(:,2))]);
        ylabel('Trials')
        title(vnames{pm})
    end
end



% 5. How does performance develop across stages?
% a. Performance parameters across all trials
nf = 5;
figure(nf+1)
comptype = 'prwin'; % Define which parameter
dz = cell(1,1); % Initiate array
lspec = {'-', ':'}; % Line style in plot

for gr = 3
    for st = 1:9 %unique(behavior.singleTrials{gr}{ra,trn}(:,3))'
        for ra = 1:size(behavior.accwin{gr},1)
            subplot(2,1,1)
            stcount = sum(behavior.singleTrials{gr}{ra,1}(:,3)==st);
            plot(st-.5:1/(stcount-1):st+.5, behavior.(comptype){gr}{ra,1}(behavior.singleTrials{gr}{ra,1}(:,3)==st)+ra*1.33,...
                'color', col_rep(st),'linewidth',2, 'linestyle', lspec{1})
            hold on
            if ra == 1
               ylabel('Rats') 
               xlabel('Stages')
            end
            set(gca,'box','off')
            
            % categorize into bins
            subplot(2,1,2)
            bincount = 20;
            dt = st-.5:1/(stcount-1):st+.5;
            dy = behavior.(comptype){gr}{ra,1}(behavior.singleTrials{gr}{ra,1}(:,3)==st);
            [~,~,binid]=histcounts(dt,bincount);
            if isempty(binid)
                dz{gr}{st,1}(ra,:) = nan(1,bincount);
            else
                for z = 1:max(binid)
                    dz{gr}{st,1}(ra,z) =  mean(dy(binid==z));
                end
            end
            
        end
        dz{gr}{st,1}(dz{gr}{st,1}==0)=nan;
        jbfill(st-.5:1/(bincount-1):st+.5,...
            nanmean(dz{gr}{st,1})+nanstd(dz{gr}{st,1})/sqrt(size(dz{gr}{st,1},1)),...
            nanmean(dz{gr}{st,1})-nanstd(dz{gr}{st,1})/sqrt(size(dz{gr}{st,1},1)),...
            col_rep(st), col_rep(st));
        hold on
        plot(st-.5:1/(bincount-1):st+.5, nanmean(dz{gr}{st,1}),...
            'color', col_rep(st),'linewidth',2, 'linestyle', lspec{1}) %Means
        hold on
        if st == 9
            xlabel('Stages')
            ylabel('% (behavioral parameter)')
        end
    end
end


% b. Distinction between trained and untrained trials
nf = 8;
figure(nf+1)
comptype = 'accwin'; % Choose from 'accwin' (accuracy), 'omwin (omissions), or 'prwin' (prematures)
dz = cell(1,1);
lspec = {'-', ':'};
bincount = 30; % Desired bin count

for gr = 3 % Loop through groups (1=only dmPFC, 2=only vmPFC, 3=all)
    for trn=1:2 % Loop through training status (untrained / trained)
        for st = 1:9 %unique(behavior.singleTrials{gr}{ra,trn}(:,3))' % Loop through stages (1= stage 1, 11= vITI)
            for ra = 1:size(behavior.accwin{gr},1) % Loop through rats
                subplot(2,1,1)
                stcount = sum(behavior.singleTrials{gr}{ra,trn+1}(:,3)==st); % trialcount
                plot(st+.5*(trn-1):.5/(stcount-1):st+.5*trn, behavior.(comptype){gr}{ra,trn+1}(behavior.singleTrials{gr}{ra,trn+1}(:,3)==st)+ra*1.33,...
                    'color', col_rep(st),'linewidth',2, 'linestyle', lspec{trn}) % Plot individual rat sliding performance window
                hold on
                if ra == 1
                    ylabel('Rats')
                    xlabel('Stages')
                end
                
                % categorize into bins
                subplot(2,1,2)
                dt = st:1*trn/(stcount-1):st+1; % x-axis values (so that trained and untrained sections each comprise 0.5 units on x-axis
                dy = behavior.(comptype){gr}{ra,trn+1}(behavior.singleTrials{gr}{ra,trn+1}(:,3)==st); % Behavioral parameter values
                [~,~,binid]=histcounts(dt,bincount); % Divide trials into bins
                if isempty(binid)
                    dz{gr}{st,trn+1}(ra,:) = nan(1,bincount); % NaN if there are no trials in this category
                else
                    for z = 1:max(binid) % Loop through bins
                        dz{gr}{st,trn+1}(ra,z) =  mean(dy(binid==z)); % Store mean traces 
                    end
                end
            end
        dz{gr}{st,trn+1}(dz{gr}{st,trn+1}==0)=nan;
        jbfill(st+.5*(trn-1)-1:.5/(bincount-1):st+.5*trn-1,...
            nanmean(dz{gr}{st,trn+1})+nanstd(dz{gr}{st,trn+1})/sqrt(size(dz{gr}{st,trn+1},1)),...
            nanmean(dz{gr}{st,trn+1})-nanstd(dz{gr}{st,trn+1})/sqrt(size(dz{gr}{st,trn+1},1)),...
            col_rep(st), col_rep(st)); % SEM confidence intervals
        hold on
        plot(st+.5*(trn-1)-1:.5/(bincount-1):st+.5*trn-1, nanmean(dz{gr}{st,trn+1}),...
            'color', col_rep(st),'linewidth',2, 'linestyle', lspec{trn}) % Plot mean traces
        hold on
        xlim([0 9])
        if st == 11 && trn ==2
            xlabel('Stages')
            ylabel('% (behavioral parameter)')
            yl = get(gca,'ylim');
            line([1 2], [yl(2)* 1.1 yl(2)*1.1], 'color', 'k', 'linewidth',2)
            line([1 2], [yl(2)*1.15 yl(2)*1.15], 'color', 'k', 'linewidth',2, 'linestyle', ':')
            text(2.1, yl(2)*1.1, 'Untrained')
            text(2.1, yl(2)*1.15, 'Trained')
            set(gca, 'ylim', [yl(1), yl(2)*1.2])
        end
        end
    end
end

% c. Difference between trained/untrained within and across stages
% pst_tr = cell(1,1);
nf = 9;
figure(nf+1)
cst_utr = cell(1,1);
cst_tr = cell(1,1);
ctypes = {'accwin' 'omwin' 'prwin'};
plotnames = {'Accuracy' 'Omissions' 'Prematures'};
lims = [1.05 0.6 0.4];
mclines = [-.17 -.03; 0.03 .17; -.17 .17]; %x-coords for multiple comparison significance lines
mclims = [1; 1; 1.05]; %y-coords for mc significance lines
utr_means = cell(1,1);
tr_means = cell(1,1);
test_results = cell(1,1);

for gr = 3 % Loop through groups (1=only dmPFC, 2=only vmPFC, 3=all)
    for ct = 1:numel(ctypes)
        for st = 1:9 %unique(behavior.singleTrials{gr}{ra,trn}(:,3))' % Loop through stages (1= stage 1, 11= vITI)
            for ra = 1:size(behavior.accwin{gr},1) % Loop through rats
                %             pst_utr{ra,st} = behavior.(comptype){gr}{ra,2}(behavior.singleTrials{gr}{ra,2}(:,3)==st-1); % Previous stage trained trials
                %             pst_tr{ra,st} = behavior.(comptype){gr}{ra,3}(behavior.singleTrials{gr}{ra,3}(:,3)==st-1); % Previous stage trained trials
                cst_utr{ra,st} = behavior.(ctypes{ct}){gr}{ra,2}(behavior.singleTrials{gr}{ra,2}(:,3)==st); % Current stage untrained trials
                cst_tr{ra,st} = behavior.(ctypes{ct}){gr}{ra,3}(behavior.singleTrials{gr}{ra,3}(:,3)==st); % Current stage trained trials
                
                % Store mean rat values
                try
                    utr_means{ct}(ra,st) = nanmean(cst_utr{ra,st});
                    tr_means{ct}(ra,st) = nanmean(cst_tr{ra,st});
                catch
                    utr_means{ct}(ra,st) = nan(1,1);
                    tr_means{ct}(ra,st) = nan(1,1);
                end
                st_diff(ra,st) = nanmean(cst_tr{ra,st})-nanmean(cst_utr{ra,st});
                
            end
            
            if st==3 || st==9
                perf_tab = [tr_means{ct}(:,st-1), utr_means{ct}(:,st) tr_means{ct}(:,st)]; % Here you can choose the comparison
                % utr_means = untrained means
                % tr_means = trained means
                
                perf_tab(sum(isnan(perf_tab),2)>0,:)=[];
                subplot(3,1,ct)
                plot(ones(size(perf_tab,1),1)*st-.2, perf_tab(:,1), '.', 'color', 'k', 'markersize', 12)
                hold on
                plot(ones(size(perf_tab,1),1)*st, perf_tab(:,2), '.', 'color', 'b', 'markersize', 12)
                plot(ones(size(perf_tab,1),1)*st+.2, perf_tab(:,3), '.', 'color', 'k', 'markersize', 12)

                plot(repmat((.8:.2:1.2)+st-1,size(perf_tab,1),1)', perf_tab',...
                    'color','k')
                
                if size(perf_tab,2)==2
                    [~,p,~,test_out]=ttest(perf_tab(:,1), perf_tab(:,2));
                    test_results{ct}(st,1:3)=[test_out.tstat, test_out.df, p];
                else
                    [p,Tab,test_out]=anova1(perf_tab,[],'off');
                    test_results{ct}(st,1:3)=[Tab{2,5}, Tab{2,3}, Tab{2,6}];
                end
                yl = get(gca, 'ylim');
                if p<0.05
                    if size(perf_tab,2)>2
                       mc=multcompare(test_out,[],'off'); 
                       test_results{ct}(st,4:6)=mc(:,6)';
                       plot(st+mclines(mc(:,6)<0.05,:)', repmat(lims(ct).*mclims(mc(:,6)<0.05,:)',2,1), 'color','k', 'linewidth',2)
%                        text(st, lims(ct), '*', 'fontsize', 14, 'horizontalalignment','center')
                    else
                        text(st, lims(ct), '*', 'fontsize', 14, 'horizontalalignment','center')
                    end
                end
                if st == 9
                    set(gca, 'ylim', [yl(1), lims(ct)+0.1*lims(ct)])
                    set(gca, 'box','off')
                end
            end
        end
        title(plotnames{ct})
        ylabel('Proportion')
        xlabel('Stage')
    end
    supertitle(['Performance difference between stages (tr-1:tr)' newline])
end


% 6. vSD and vITI
nf = nf+2;
vnames = {'accuracy' 'omissionpct' 'prematurepct'};
figure(nf);
statmat = cell(3,2);
test_results = cell(3,2);

for pm = 1:numel(vnames)
    subplot(numel(vnames),1,pm)
    
    % dmPFC (s9) - dmPFC(easySD midSD hardSD)
    
    behavior.(vnames{pm}){1}(behavior.(vnames{pm}){1}==0)=nan;
    behavior.(vnames{pm}){2}(behavior.(vnames{pm}){2}==0)=nan;
    behavior.(vnames{pm}){3}(behavior.(vnames{pm}){3}==0)=nan;
    hold on
    
    for gr = 2:3
        
        plot(ones(numel(behavior.(vnames{pm}){gr}(9,:,1)))*gr, behavior.(vnames{pm}){gr}(9,:,1), '.', 'color','k', 'markersize', 12)
        statmat{pm,gr}(1,:) = behavior.(vnames{pm}){gr}(9,:,1);
        for ss = 1:3
            plot(ones(numel(behavior.(vnames{pm}){gr}(10,:,ss)))*gr+.125*ss, behavior.(vnames{pm}){gr}(10,:,ss), '.', 'color',col_rep(gr), 'markersize', 12)
            statmat{pm,gr}(ss+1,:) = behavior.(vnames{pm}){gr}(10,:,ss);
        end
        
        % stats
        [P,TAB,STR]=anova1(statmat{pm,gr}', [], 'off');
        
        if P < 0.05
            [MC1, ~] = multcompare(STR, 'display', 'off');
            test_results{pm,gr}=[P; TAB{2,3}; TAB{2,5}; MC1(:,end)];
            text(gr-0.1,  nanmean(behavior.(vnames{pm}){gr}(9,:,1)), '***','HorizontalAlignment', 'center')
        end
                 
    end
    
    ylabel('Proportion')
    title(vnames{pm})
    set(gca,'box','off')
    
    xlim([0.8 2.5])
    
end

% Save figures
figHandles = findobj('Type', 'figure'); % Gather all figure handles
set(figHandles, 'renderer', 'painters'); % Set renderer to 'painters' for proper further processing in Illustrator
set(figHandles, 'Position', get(0, 'Screensize')); % Set figures to full screen size (easier to process images in Illustrator)

for fi = 1:numel(figHandles)
    saveas(figHandles(fi), [savefolder2, '\Behavior-',num2str(figHandles(fi).Number)], 'epsc'); % Save each figure as .EPSC file 
end


%% Figures (vITI and vSD included)
close all

% 6. Behavioral performance in variable ITI and SD sessions compared to
% baseline (stage 9 trained trials)
nf = 9;
figure(nf+1)
vnames9 = {'accuracy' 'omissionpct' 'prematurepct' 'correctlatency' 'incorrectlatency' 'posterrorlatency' 'prematuretime'};
gr = 3;

for vn = 1:numel(vnames9)
    
    % Upper row of plots shows Stage 9 and vSD and vITI with different ITI
    % and SD conditions
    subplot(3,numel(vnames9),vn)
    hold on
    % Stage 9 boxplot
    atmp = behavior.(vnames9{vn}){gr}(9,:,1);
    atmp(atmp==0)=nan;
    plot(ones(size(behavior.(vnames9{vn}){gr}(9,:,1),2),1)', atmp,'.','markersize', 14, 'color',col_resp(vn))
%     boxplot(behavior.(vnames9{vn}){gr}(9,:,1), 'positions', 1,'plotstyle','compact','colors',col_resp(vn))
    yl(1,:) = get(gca,'ylim');
    
    %vSD boxplot
    btmp=reshape(behavior.(vnames9{vn}){gr}(10,:,1:3), [size(behavior.(vnames9{vn}){gr}(10,:,1:3),2),3]);
    btmp(btmp==0)=nan;
    plot(repmat((1.75:0.25:2.25)',1, size(behavior.(vnames9{vn}){gr}(10,:,1:3),2)),btmp', '.','markersize', 10, 'color',col_resp(vn)) % data points
    plot(repmat((1.75:0.25:2.25)',1, size(behavior.(vnames9{vn}){gr}(10,:,1:3),2)),btmp', 'color',col_resp(vn), 'linewidth', 1) %Lines
    yl(2,:) = get(gca,'ylim');
    
    %vITI boxplot
    ctmp=reshape(behavior.(vnames9{vn}){gr}(11,:,1:3), [size(behavior.(vnames9{vn}){gr}(11,:,1:3),2),3]);
    ctmp(ctmp==0)=nan;
    plot(repmat((2.75:0.25:3.25)',1, size(behavior.(vnames9{vn}){gr}(10,:,1:3),2)),ctmp', '.','markersize', 10, 'color',col_resp(vn)) % Data points
    plot(repmat((2.75:0.25:3.25)',1, size(behavior.(vnames9{vn}){gr}(10,:,1:3),2)),ctmp', 'color',col_resp(vn), 'linewidth', 0.5) %Lines
    yl(3,:) = get(gca,'ylim');
    
    % Tests
    [vsd_p,~,vsd_stats]=anova1(btmp,[],'off');
    if vsd_p<0.05
        vsd_mc = multcompare(vsd_stats,[],'off');
        text(2, max(yl(:,2)), '*', 'fontsize', 18, 'fontweight', 'b')
    end
    
    [viti_p,~,viti_stats]=anova1(ctmp,[],'off');
    if viti_p<0.05
        viti_mc = multcompare(viti_stats,[],'off');
        text(3, max(yl(:,2)), '*', 'fontsize', 18, 'fontweight', 'b')
    end
    
    % Plot parameters
    xlim([0.5 3.5])
    title(vnames9{vn})
    ylabel('Proportion')
    newyl = [max([0,min(yl(:,1))]), max(yl(:,2))+0.1*max(yl(:,2))];
    line([1.5 1.5],newyl, 'color','k','linestyle','-.')
    line([2.5 2.5],newyl, 'color','k','linestyle','-.')
    set(gca,'box','off')
    xlabel('Session Type')
    set(gca,'xtick', 1:3, 'xticklabels', {'s9' 'vSD' 'vITI'}, 'XTickLabelRotation',45)
    set(gca, 'ylim', newyl)
    
    % Middle row compares stage 9 with average of all task parameters in vSD and vITI
    % sessions
    
    subplot(3,numel(vnames9),numel(vnames9)+vn)
    hold on
    all_trials = [behavior.(vnames9{vn}){3}(9,:,1);...
        nanmean(behavior.(vnames9{vn}){gr}(10,:,:),3);...
        nanmean(behavior.(vnames9{vn}){gr}(11,:,:),3)]; % Matrix with stage 9, vSD(1s cue), and vITI(5s delay)
    all_trials(all_trials==0)=nan;

    [bt_p,~,bt_stats]=anova1(all_trials',[],'off');
    plot(repmat((1:3),size(all_trials,2),1)',all_trials, '.', 'markersize', 10,'color', col_resp(vn)); % Markers
    plot(repmat((1:3),size(all_trials,2),1)',all_trials, 'linewidth', 0.5,'color', col_resp(vn)); % Line
    
    yl2 = get(gca, 'ylim');
    if bt_p<.05
        bt_mc = multcompare(bt_stats,[],'off');
        text(2, yl2(2)+0.1*yl2(2), '*', 'fontsize', 18, 'fontweight', 'b')
    end
    set(gca, 'ylim', [yl2(1) yl2(2)*0.15+yl2(2)])
    xlim([0.5 3.5])
    ylabel('Proportion')
    xlabel('Session Type')
    set(gca,'xtick', 1:3, 'xticklabels', {'s9' 'vSD' 'vITI'}, 'XTickLabelRotation',45)
    
    % Lower row compares stage 9 with same task parameters in vSD and vITI
    % sessions (1s cue, 5s delay)
    
    subplot(3,numel(vnames9),numel(vnames9)*2+vn)
    hold on
    baseline_trials = [behavior.(vnames9{vn}){3}(9,:,1);...
        behavior.(vnames9{vn}){gr}(10,:,1);...
        behavior.(vnames9{vn}){gr}(11,:,1)]; % Matrix with stage 9, vSD(1s cue), and vITI(5s delay)
    baseline_trials(baseline_trials==0)=nan;

    [bt_p,~,bt_stats]=anova1(baseline_trials',[],'off');
    plot(repmat((1:3),size(baseline_trials,2),1)',baseline_trials, '.', 'markersize', 10,'color', col_resp(vn));
    plot(repmat((1:3),size(baseline_trials,2),1)',baseline_trials, 'linewidth', 0.5, 'color', col_resp(vn));
       yl3 = get(gca, 'ylim');
    
    if bt_p<.05
       bt_mc = multcompare(bt_stats,[],'off');
       text(2, yl3(2)+0.1*yl3(2), '*', 'fontsize', 18, 'fontweight', 'b')
    end
    set(gca, 'ylim', [yl3(1) yl3(2)*0.15+yl3(2)])
    xlim([0.5 3.5])
    ylabel('Proportion')
    xlabel('Session Type')
    set(gca,'xtick', 1:3, 'xticklabels', {'s9' 'vSD' 'vITI'}, 'XTickLabelRotation',45)
end


% 7. Progression of performance in vITI/vSD sessions. 
% 



% N. Predictive power of learning stages in sessions that require higher
% cognitive load.
%
% Divide into good learners and bad learners? 
%   - least trials needed
%       + Highest vs lowest half of rats?
%       1) Total number of trials
%       2) Rating per stage, see which stage is most predictive
%   - final performance in stage 9  
%       + need to selectively look at last N trials (N=50? 100?), to ensure
%       final baseline performance
% Make plots for e.g. correct responses, different SDs
vnames = {'accuracy', 'omissionpct', 'prematurepct', 'correctlatency', 'prematuretime', 'posterrorlatency'};
msize = 6; 
nf = 10; % Number of figures

for vn = 1:numel(vnames)
    tmp1 = [behavior.(vnames{vn}){1}(10,:,1); behavior.(vnames{vn}){1}(10,:,2); behavior.(vnames{vn}){1}(10,:,3)]; % vSD performance (accuracy)
    tmp2 = [behavior.(vnames{vn}){1}(11,:,1); behavior.(vnames{vn}){1}(11,:,2); behavior.(vnames{vn}){1}(11,:,3)]; % vITI performance (accuracy)
    tmp1null = isnan(tmp1(1,:));
    tmp2null = isnan(tmp2(1,:));
    % highCogSD = (nanmean(tmp1)>nanmedian(nanmean(tmp1),2))'; % High and low performing rats in vSD
    % highCogITI = (nanmean(tmp2)>nanmedian(nanmean(tmp2),2))'; % High and low performing rats in vITI
    
    
    % Average of vSD/ITI performance for each rat
    tmp1=tmp1(:,~tmp1null);
    tmp2=tmp2(:,~tmp2null);
    
    % Per stage, relation between performance in vITI and vSD subconditions and learning
    % speed (number of trials required to finish each stage)
    figure(nf+1)
    for stage = 1:9
        subplot(numel(vnames),2,2*vn-1) % var SD
        learnRatingSD = behavior.trialsBeforeTrained{1}(~tmp1null,:)<=nanmedian(behavior.trialsBeforeTrained{1}(~tmp1null,:)); % Number of trials required per stage
        plot(stage+(-0.25:0.2:0.15).*ones(sum(learnRatingSD(:,stage)==1),1), tmp1(:,learnRatingSD(:,stage)==1)', '.b','markersize',msize)
        hold on
        bar(stage+(-0.28:0.2:0.12),mean(tmp1(:,learnRatingSD(:,stage)==1),2), 0.15, 'facecolor','b', 'edgecolor', 'b',  'facealpha', 0.4); % bar with means
        plot(stage+(-0.15:0.2:0.25).*ones(sum(learnRatingSD(:,stage)==0),1), tmp1(:,learnRatingSD(:,stage)==0)', '.r','markersize',msize)
        bar(stage+(-0.19:0.2:0.21),mean(tmp1(:,learnRatingSD(:,stage)==0),2), 0.15, 'facecolor','r', 'edgecolor', 'r',  'facealpha', 0.4); % bar with means
        
        if stage == 1
            ylabel(vnames{vn})
        end
        if vn == 1
            title('vSD')
        end
        
        learnRatingITI = behavior.trialsBeforeTrained{1}(~tmp2null,:)<=nanmedian(behavior.trialsBeforeTrained{1}(~tmp2null,:)); % Number of trials required per stage
        subplot(numel(vnames),2,2*vn) % Var ITI
        plot(stage+(-0.25:0.2:0.15).*ones(sum(learnRatingITI(:,stage)==1),1), tmp2(:,learnRatingITI(:,stage)==1)', '.b','markersize',msize)
        hold on
        bar(stage+(-0.28:0.2:0.12),mean(tmp2(:,learnRatingITI(:,stage)==1),2), 0.15, 'facecolor','b', 'edgecolor', 'b', 'facealpha', 0.4); % bar with means
        plot(stage+(-0.15:0.2:0.25).*ones(sum(learnRatingITI(:,stage)==0),1), tmp2(:,learnRatingITI(:,stage)==0)', '.r','markersize',msize)
        bar(stage+(-0.19:0.2:0.21),mean(tmp2(:,learnRatingITI(:,stage)==0),2), 0.15, 'facecolor','r', 'edgecolor', 'r',  'facealpha', 0.4); % bar with means
        
        if stage == 1
            ylabel(vnames{vn})
        end
        if vn == 1
            title('vITI')
        end
    end
    
    
    % Split by subcondition (SDs/ ITIs)
    figure(nf+2)
    for stage = 1:9
        learnRatingSD = behavior.trialsBeforeTrained{1}(~tmp1null,:)>nanmedian(behavior.trialsBeforeTrained{1}(~tmp1null,:)); % Number of trials required per stage
        % Signs:
        % '>' : BLUE lines mean slower learners (need more trials than median
        % of all rats in stage)
        subplot(2,numel(vnames),vn) % var SD
        sdtmp_1(stage,:) = mean(tmp1(:,learnRatingSD(:,stage)==1),2);
        sdtmp_2(stage,:) = mean(tmp1(:,learnRatingSD(:,stage)==0),2);
        plot(sdtmp_1(stage,:),'b')
        hold on
        plot(sdtmp_2(stage,:),'r')
        if stage == 9
            title(vnames{vn})
        end
        
        learnRatingITI = behavior.trialsBeforeTrained{1}(~tmp2null,:)>nanmedian(behavior.trialsBeforeTrained{1}(~tmp2null,:)); % Number of trials required per stage
        subplot(2,numel(vnames),numel(vnames)+vn) % Var ITI
        ititmp_1(stage,:)=mean(tmp2(:,learnRatingITI(:,stage)==1),2);
        ititmp_2(stage,:)=mean(tmp2(:,learnRatingITI(:,stage)==0),2);
        plot(ititmp_1(stage,:),'b')
        hold on
        plot(ititmp_2(stage,:),'r')
    end
end

% Save figures
figHandles = findobj('Type', 'figure'); % Gather all figure handles
set(figHandles, 'renderer', 'painters'); % Set renderer to 'painters' for proper further processing in Illustrator
set(figHandles, 'Position', get(0, 'Screensize')); % Set figures to full screen size (easier to process images in Illustrator)

for fi = 1:numel(figHandles)
    saveas(figHandles(fi), [savefolder2, '\Behavior-',num2str(figHandles(fi).Number)], 'epsc'); % Save each figure as .EPSC file 
end


%% 3) Signal analysis
% This chapter will analyze the calcium signals corresponding to distinct
% 5-CSRTT stages.
%
% Still under construction!
%
% Input: 
%   a) meta matrix called 'meta_in', generated in chapter 1
%   b) data matrix called 'data_in(_z)', generated in chapter 1
%
% Output: 
%
%   3.0) Structure called 'calcium', which contains traces and signal
%   parameters
%
% Figures of different behavioral analyses:
%
% 3.11) Means of rats, per stage, separate figures
% 3.11a) Means of rats, stage 3 vs stage 9, with permutation test
% 3.11b) Comparisons between stage 3 (early, putatively unshaped activity)
% and stage 9 (late, shaped activity)
% 3.11c) Signal parameters per stage, comparisons
% Compares signal parameters across stages
% 3.12) Means of rats, per stage, one figure per group. 
% 3.12a) Means of rats, stage by stage comparison  
% 3.12b) Means of rats, correlations between stages
% 3.13) Means of rats, per stage, separate figures, outcomes split
% 3.14) Means of rats, per stage, one figure per group, split by outcomes
% 3.15) Means of rats, per stage, separate figures, outcomes AND trained
% status split
% 3.16a) Means of rats, per stage, separate figures, dmPFC vs vmPFC
% 3.16b) Signal parameters, dmPFC vs vmPFC
% 3.17a) How does signal develop across stages? Similar figure to behavior
% figure 2.6a
% 3.17b). Distinction between trained and untrained trials. Similar to
% 3.17a, but now trained and untrained trials are separated

%
% Other useful things:
% Meta table contents (useful when tweaking this chapter)
% Generate or append trial data matrix
% 1      2    3     4     5     6        7           8
% group, rat, cond, sess, resp, subspec, lower norm, upper norm
% 9             % 10         % 11
% trial time,   resp time    learning/nolearning

%% 3.0) Generate traces and signal parameters
% Initiate variables
win{1} = ceil(5*15.89):ceil(15*15.89); % ITI window + cue + first responses
win{2} = ceil(6.5*15.89):ceil(16.5*15.89); % Response window
bhvwin{1} = ceil(5*15.89): ceil(10*15.89); % ITI window
bhvwin{2} = ceil(15.89):ceil(5*15.89); % Around response
calcium = struct;

for sy = 1:2 % Loop through sync moments
    % Initiate signal storage structure
%     data = (data_in{sy}(:,win{sy})-meta_in(:,7))./(meta_in(:,8)-meta_in(:,7)); % Proportionalized data to 1st and 99th percentile
    data = data_in_z{sy}; % Non-proportionalized
    meta = meta_in;
    
    % Calculate AUC, information distance, other parameters
    calcium(sy).traces = data; % Regular calcium traces (proportionalized between 1-99 percentile)
    calcium(sy).auc = trapz(data(:,bhvwin{sy}),2); %AUC
    [calcium(sy).maxA, calcium(sy).maxT] = max(data(:,bhvwin{sy}),[],2); % Maximum and timing
    calcium(sy).amplitude = max(data(:,bhvwin{sy}),[],2)-min(data(:,bhvwin{sy}),[],2); % Amplitude during response window
    calcium(sy).variance = var(data,[],2); % Variance
    if sy==1
        calcium(1).dfiti = data(:,win{1}(1))-data(:,win{1}(ceil(5*15.89))); % difference between sig at start delay and cue onset
        calcium(1).fcue = nanmean(data(:,win{1}(ceil(5*15.89))-5:win{1}(ceil(5*15.89))),2); % calcium signal at cue onset (compared to baseline)
        calcium(1).fstart = nanmean(data(:,(ceil(5*15.89)):(ceil(5*15.89))+5),2); % calcium signal at trial start (compared to baseline)
    elseif sy == 2
        calcium(2).fresp = data(:,win{2}(ceil(3*15.89))); % calcium signal at response (compared to baseline)
        calcium(2).peakrespdf = max(data(:,ceil(3*15.89):ceil(6*15.89)),[],2)-data(:,win{2}(ceil(3*15.89))); % difference between peak sig and pre response
        calcium(2).respdecaydf = data(:,win{2}(ceil(3*15.89)))-min(data(:,ceil(3*15.89):ceil(6*15.89)),[],2); % difference between response onset and trough after peak
        calcium(2).preresp = nanmean(data(:,ceil(15.89):ceil(3*15.89)),2); % Mean values 2s before response
    end
    
    % Specific peak values
%     calcium(sy).peakwidth = nan(size(data,1),20);
%     calcium(sy).peaktime = nan(size(data,1),20);
%     calcium(sy).peakprominence = nan(size(data,1),20);
%     calcium(sy).peakauc = nan(size(data,1),20);
%     for tr = 1:size(data_in{1},1)
%         [PK, LO, WI, PR] = findpeaks(data(tr,bhvwin{sy}), 'MinPeakProminence', 0.2*(max(data(tr,bhvwin{sy}))-min(data(tr,bhvwin{sy})))); % Get peaks that are at least 20% of max amplitude
%         calcium(sy).peakwidth(tr,1:numel(WI)) = WI;% Peak width
%         calcium(sy).peaktime(tr,1:numel(LO)) = LO;% Peak timing
%         calcium(sy).peakprominence(tr,1:numel(PR)) = PR;% Peak prominence
%         calcium(sy).npeaks(tr,1) = numel(PK); % Number of peaks
%         b_min=[];b_max=[];
%         for k = 1:length(LO)% Peak AUC
%             i = LO(k);
%             while i > 1 && data(tr,i-1) > (PK(k)-PR(k))
%                 i = i - 1;
%             end
%             b_min(k) = i;
%             i = LO(k);
%             while i < length(data(tr,:)) && data(tr,i+1) > (PK(k)-PR(k))
%                 i = i + 1;
%             end
%             b_max(k) = i;
%             calcium(sy).peakauc(tr,k) = trapz(data(tr,b_min(k):b_max(k))-(PK(k)-PR(k)));
%         end
%     end
    
%     % Aggregated specific peak values
%     calcium(sy).meanpeakwidth = nanmean(calcium(sy).peakwidth,2);
%     calcium(sy).meanpeakprominence = nanmean(calcium(sy).peakprominence,2);
%     calcium(sy).meanpeakauc = nanmean(calcium(sy).peakauc,2);
%     calcium(sy).totalpeakauc = nansum(calcium(sy).peakauc,2);
%     
    % Calcium raw activity during iti window
    calcium(sy).weight = trapz(data(:,bhvwin{sy})-min(data(:,bhvwin{sy}),[],2),2); % Mean bulk activity (measured from 1st percentile of signal in entire session)
    calcium(sy).mean = nanmean(data(:,bhvwin{sy}),2); % Mean values (based on trial-related baseline)
    calcium(sy).slope = (nanmean(data(:,bhvwin{sy}(1:5)),2)-nanmean(data(:,bhvwin{sy}(end-5:end)),2))/5; % Slope during delay (same as dfiti?)
    
    if sy == 1
        %Raw
        calcium(sy).weight1 = trapz(data(:,bhvwin{sy}(1:round(.5*length(bhvwin{sy}))))-meta_in(:,7),2); %min(data(:,bhvwin{sy}),[],2),2); % During first half of ITI
        calcium(sy).weight2 = trapz(data(:,bhvwin{sy}(round(.5*length(bhvwin{sy})):end))-meta_in(:,7),2); %min(data(:,bhvwin{sy}),[],2),2); % During second half of ITI
        calcium(sy).mean1 = nanmean(data(:,bhvwin{sy}(1:round(.5*length(bhvwin{sy})))),2); % Mean values during first half of ITI
        calcium(sy).mean2 = nanmean(data(:,bhvwin{sy}(round(.5*length(bhvwin{sy})):end)),2); % Mean values during second half of ITI
        calcium(sy).slope1 = nanmean(data(:,bhvwin{sy}(1:5)),2)-nanmean(data(:,bhvwin{sy}(round(end*0.5)-5:round(end*0.5))),2); % Slope during delay (same as dfiti?)
        calcium(sy).slope2 = nanmean(data(:,bhvwin{sy}(round(end*0.5)-5:round(end*0.5))),2)-nanmean(data(:,bhvwin{sy}(end-5:end)),2); % Slope during delay (same as dfiti?)
        calcium(sy).bias = calcium(sy).mean1-calcium(sy).mean2;

        %Proportion
        calcium(sy).weight1pct = calcium(sy).weight1./...
            calcium(sy).weight; % During first half of ITI
        calcium(sy).weight2pct = calcium(sy).weight2./...
            calcium(sy).weight;    
    end
    % Calculate at which point signal has reached 50% total activity during
    % iti
    [~,calcium(sy).halfweighttime]=max(cumtrapz(data(:,bhvwin{sy})-min(data(:,bhvwin{sy}),[],2),2)>0.5*calcium(sy).weight,[],2);
end


%% Figures
% Meta table contents (useful when tweaking this chapter)
% Generate or append trial data matrix
% 1      2    3     4     5     6        7           8
% group, rat, cond, sess, resp, subspec, lower norm, upper norm
% 9             % 10         % 11
% trial time,   resp time    learning/nolearning

% Initiate variables
vnames = fieldnames(calcium(sy));
gnames = {'dmPFC' 'vmPFC'};
snames = {'Sync at trialstart', 'Sync at response'};
stnames = {'Stage 1' 'Stage 2' 'Stage 3' 'Stage 4' 'Stage 5' 'Stage 6' 'Stage 7' 'Stage 8' 'Stage 9'};
xv{1} = 0:17.5/(size(calcium(1).traces,2)-1):17.5; % x-axis for combined sync plots (trialstart)
xv{2} = 18.5+(0:17.5/(size(calcium(2).traces,2)-1):17.5);  % x-axis for combined sync plots (response)

% 3.10) Heatmap of all trials
for gr = 1:1
    figure
    for sy = 1:1
        k = 1;
%         for ra = unique(meta_in(meta_in(:,1)==gr,2))'
        for ra = 10
            trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)<10 & meta_in(:,5)==1);
            imagesc(calcium(sy).traces(trid,:))
            caxis([-5 5])
            k = k+1;
        end
    end
end

% 3.11) Means of rats, per stage, separate figures
nf = 9; % number of figures in behavioral part
data_tmp = cell(1,1);

for gr = 1:2
    figure(nf+gr-1)
    for st = 1:9
        for sy = 1:2
            k=1;
            for ra = unique(meta_in(meta_in(:,1)==gr,2))'
                trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st & meta_in(:,5)==1);
                data_tmp{sy}(k,:)=nanmean(calcium(sy).traces(trid,:));
                k=k+1;
            end
            subplot(3,3,st)
            plot(xv{sy},nanmean(data_tmp{sy}), 'color', col_rep(gr), 'linewidth',2) % Plot means
            hold on
            jbfill(xv{sy}, nanmean(data_tmp{sy})+nanstd(data_tmp{sy})/sqrt(size(data_tmp{sy},1)),...
                nanmean(data_tmp{sy})-nanstd(data_tmp{sy})/sqrt(size(data_tmp{sy},1)),...
                col_rep(gr), col_rep(gr)); % Plot SEM
            hold on
            title(['stage ' num2str(st)])
            if sy == 2
                yl = get(gca, 'ylim'); % Get y-axis limits
            line([5 5], yl, 'color','k')
            line([10 10], yl, 'color','k')
            line([24.5 24.5], yl, 'color','k')
            set(gca,'ylim',yl)
            end
        end
    end
    supertitle(gnames{gr})
end

% 3.11a) Means of rats, stage 3 vs stage 9, with permutation test
nf=nf+2; %Figure indices
figure(nf)
st1 = 3;
st2 = 9;

for gr = 1:2
    subplot(2,1,gr)
    data_tmp0 = cell(1,1);
    for sy = 1:2
        k=1;
        for ra = unique(meta_in(meta_in(:,1)==gr,2))'
            trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,5)==1 & meta_in(:,3)==st1); % Only stage 3 trials
            trid2 = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,5)==1 & meta_in(:,3)==st2); % Only stage 9 trials
            data_tmp0{sy,1}(k,:)=nanmean(calcium(sy).traces(trid,:)); % Stage 3
            data_tmp0{sy,2}(k,:)=nanmean(calcium(sy).traces(trid2,:)); % Stage 9
            k=k+1;
        end
        
        % Plot graphs
        q(1) = plot(xv{sy},nanmean(data_tmp0{sy,1}), 'color', 'k', 'linewidth',2); % Plot s3 means
        hold on
        q(2) = plot(xv{sy},nanmean(data_tmp0{sy,2}), 'color', col_rep(gr), 'linewidth',2); % Plot s9 means
        jbfill(xv{sy}, nanmean(data_tmp0{sy,1})+nanstd(data_tmp0{sy,1})/sqrt(size(data_tmp0{sy,1},1)),...
            nanmean(data_tmp0{sy,1})-nanstd(data_tmp0{sy,1})/sqrt(size(data_tmp0{sy,1},1)),...
            'k', 'k'); % Plot SEM s3
        hold on
        jbfill(xv{sy}, nanmean(data_tmp0{sy,2})+nanstd(data_tmp0{sy,2})/sqrt(size(data_tmp0{sy,2},1)),...
            nanmean(data_tmp0{sy,2})-nanstd(data_tmp0{sy,2})/sqrt(size(data_tmp0{sy,2},1)),...
            col_rep(gr), col_rep(gr)); % Plot SEM s9
        hold on
        
        % Permutation test
        sig_out = permTest_array(data_tmp0{sy,1}(~isnan(data_tmp0{sy,1}(:,1)),:), data_tmp0{sy,2}(~isnan(data_tmp0{sy,2}(:,1)),:),5000);
        yl = get(gca, 'ylim'); % Get y-axis limits
        xw = xv{sy}; xw(sig_out>0.05)=nan;
        plot(xw, ones(numel(xw),1)*yl(2)*0.95, 'linewidth', 3, 'color', 'k')

        % Plot parameters
        title(gnames{gr})
        ylabel('dF/F')
        xlabel('Time (s)')
        if sy == 2
            line([5 5], yl, 'color','k')
            line([10 10], yl, 'color','k')
            line([30.5 30.5], yl, 'color','k')
            set(gca,'ylim',yl)
        end
    end
    supertitle(['Stage 3 vs stage 9 comparison' newline])
    legend(q, {'Stage 3' 'Stage 9'}, 'location', 'best', 'orientation', 'horizontal')
end

% 3.11b) Comparisons between stage 3 (early, putatively unshaped activity)
% and stage 9 (late, shaped activity)
comps = {'amplitude' 'dfiti' 'fstart' 'fcue' 'slope' 'slope1' 'slope2' ... %9
    'mean' 'mean1' 'mean2' 'halfweighttime' 'bias' 'preresp'}; % Calcium signal parameters
startst = ones(numel(comps),1)*3; % Earliest stage parameter is relevant (i.e. there is no ITI in stage 1 and 2, so dF during ITI isn't applicable)
nf=nf+1; %Figure indices
data_comp = cell(1,1);
test_outputs = cell(1,1);
T_OUT_1 = cell(1,1);

for sy = 1:2 % Loop syncs
    figure(nf+sy-1)
    for gr = 1:2 % Loop groups
        for co = 1:numel(comps)
            pval = ones(1,2);
            % Compare early and late stage trial data for each parameter
            k=1;
            for ra = unique(meta_in(meta_in(:,1)==gr,2))' % Loop rats
                trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==4 & meta_in(:,5)==1); % Select trials
                trid2 = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==3 & meta_in(:,5)==1); % Select trials
                
                try % By-pass empty comparisons (sometimes it depends on sync whether or not there is one)
                    data_comp{gr}{sy,co}(k,:) = ...
                        [nanmean(calcium(sy).(comps{co})(trid,:)), nanmean(calcium(sy).(comps{co})(trid2,:))]; % Store mean values for each stage
                catch
                    data_comp{gr}{sy,co}(k,:)=nan(1,2);
                end
                k=k+1;
            end
            
            subplot(ceil(sqrt(numel(comps))), round(sqrt(numel(comps))), co) % Subplot
            title(comps{co})
            plot(repmat(1, size(data_comp{gr}{sy,co},1), 1)+3*(gr-1),data_comp{gr}{sy,co}(:,1),'.','color','k', 'markersize', 10) % Scatter plot of all values
            hold on
            plot(repmat(2, size(data_comp{gr}{sy,co},1), 1)+3*(gr-1),data_comp{gr}{sy,co}(:,2),'.','color',col_rep(gr), 'markersize', 10) % Scatter plot of all values
            plot((1:2)+3*(gr-1), nanmean(data_comp{gr}{sy,co}), 'xk') % Plot mean
            plot(repmat(1:2, size(data_comp{gr}{sy,co},1),1)'+3*(gr-1),data_comp{gr}{sy,co}', 'color','k') % Plot individual rat traces
            set(gca, 'xlim', [0 6])
            xlabel('<-- dmPFC || vmPFC -->')

            % Test
            [~,pval(gr),~,s_out]=ttest(data_comp{gr}{sy,co}(:,1), data_comp{gr}{sy,co}(:,2));
            test_outputs{gr,sy}(:,co) = [round(s_out.tstat,4), s_out.df, round(pval(gr),4)];
            if gr == 2
                yl = get(gca, 'ylim');
                for gr2 = 1:2
                    if pval(gr2)<0.05
                        line([1 2] +3*(gr2-1), [ yl(2)+(yl(2)-yl(1))*0.1  yl(2)+(yl(2)-yl(1))*0.1], 'color', 'k')
                        text(1.5+3*(gr2-1),  yl(2)+(yl(2)-yl(1))*0.1, '*','HorizontalAlignment', 'center')
                        if pval(gr2)<0.01
                            text(1.5+3*(gr2-1),  yl(2)+(yl(2)-yl(1))*0.1, '**','HorizontalAlignment', 'center')
                            if pval(gr2)<0.001
                                text(1.5+3*(gr2-1),  yl(2)+(yl(2)-yl(1))*0.1, '***','HorizontalAlignment', 'center')
                            end
                        end
                    end
                end
                set(gca, 'ylim', [yl(1), yl(2)+(yl(2)-yl(1))*0.2])
            end
        end
        clear pval
    T_OUT_1{gr,sy}=array2table(test_outputs{gr,sy}, 'VariableNames', comps);
    end
    supertitle(snames{sy})
end


% 3.11c) Means of rats, dmPFC vs vmPFC, with permutation test
nf=nf+2; %Figure indices
figure(nf)
stno = [4 3];
stnames = {'Stage 3' 'Stage 9'};

for st = 1:2
    subplot(2,1,st)
    data_tmp0 = cell(1,1);
    for sy = 1:2
        k=1;
        for ra = unique(meta_in(:,2))'
            trid = (meta_in(:,1)==1 & meta_in(:,2)==ra & meta_in(:,5)==1 & meta_in(:,3)==stno(st)); % Only stage 3 trials
            trid2 = (meta_in(:,1)==2 & meta_in(:,2)==ra & meta_in(:,5)==1 & meta_in(:,3)==stno(st)); % Only stage 9 trials
            data_tmp0{sy,1}(k,:)=nanmean(calcium(sy).traces(trid,:)); % dmPFC
            data_tmp0{sy,2}(k,:)=nanmean(calcium(sy).traces(trid2,:)); % vmPFC
            k=k+1;
        end
        
        % Plot graphs
        q(1) = plot(xv{sy},nanmean(data_tmp0{sy,1}), 'color', col_rep(1), 'linewidth',2); % Plot means
        hold on
        q(2) = plot(xv{sy},nanmean(data_tmp0{sy,2}), 'color', col_rep(2), 'linewidth',2); % Plot means
        jbfill(xv{sy}, nanmean(data_tmp0{sy,1})+nanstd(data_tmp0{sy,1})/sqrt(size(data_tmp0{sy,1},1)),...
            nanmean(data_tmp0{sy,1})-nanstd(data_tmp0{sy,1})/sqrt(size(data_tmp0{sy,1},1)),...
            col_rep(1), col_rep(1)); % Plot SEM
        hold on
        jbfill(xv{sy}, nanmean(data_tmp0{sy,2})+nanstd(data_tmp0{sy,2})/sqrt(size(data_tmp0{sy,2},1)),...
            nanmean(data_tmp0{sy,2})-nanstd(data_tmp0{sy,2})/sqrt(size(data_tmp0{sy,2},1)),...
            col_rep(2), col_rep(2)); % Plot SEM
        hold on
        
        % Permutation test
        sig_out = permTest_array(data_tmp0{sy,1}(~isnan(data_tmp0{sy,1}(:,1)),:), data_tmp0{sy,2}(~isnan(data_tmp0{sy,2}(:,1)),:),5000);
        yl = get(gca, 'ylim'); % Get y-axis limits
        xw = xv{sy}; xw(sig_out>0.05)=nan;
        plot(xw, ones(numel(xw),1)*yl(2)*0.95, 'linewidth', 3, 'color', 'k')

        % Plot parameters
        title(stnames{st})
        ylabel('dF/F')
        xlabel('Time (s)')
        if sy == 2
            line([5 5], yl, 'color','k')
            line([10 10], yl, 'color','k')
            line([21.5 21.5], yl, 'color','k')
            set(gca,'ylim',yl)
        end
    end
    supertitle(['vmPFC vs dmPFC comparison' newline])
    legend(q, {'dmPFC' 'vmPFC'}, 'location', 'best', 'orientation', 'horizontal')
end

% 3.11d) Comparisons between dmPFC and vmPFC during early (3) and late
% stages (9)
comps = {'amplitude' 'dfiti' 'fstart' 'fcue' 'slope' 'slope1' 'slope2' ... %9
    'mean' 'mean1' 'mean2' 'halfweighttime' 'bias' 'preresp'}; % Calcium signal parameters
startst = ones(numel(comps),1)*3; % Earliest stage parameter is relevant (i.e. there is no ITI in stage 1 and 2, so dF during ITI isn't applicable)
nf=nf+1; %Figure indices
data_comp = cell(1,1);
test_outputs = cell(1,1);
T_OUT_2 = cell(1,1);

for sy = 1:2 % Loop syncs
    figure(nf+sy-1)
    for st = 1:2 % Loop groups
        for co = 1:numel(comps)
            pval = ones(1,2);
            % Compare early and late stage trial data for each parameter
            k=1;
            for ra = unique(meta_in(:,2))' % Loop rats
                trid = (meta_in(:,1)==1 & meta_in(:,2)==ra & meta_in(:,3)==stno(st)); % Select trials
                trid2 = (meta_in(:,1)==2 & meta_in(:,2)==ra & meta_in(:,3)==stno(st)); % Select trials
                
                try % By-pass empty comparisons (sometimes it depends on sync whether or not there is one)
                    data_comp{st}{sy,co}(k,:) = ...
                        [nanmean(calcium(sy).(comps{co})(trid,:)), nanmean(calcium(sy).(comps{co})(trid2,:))]; % Store mean values for each stage
                catch
                    data_comp{st}{sy,co}(k,:)=nan(1,2);
                end
                k=k+1;
            end
            
            subplot(ceil(sqrt(numel(comps))), round(sqrt(numel(comps))), co) % Subplot
            title(comps{co})
            plot(repmat(1, size(data_comp{1}{sy,co},1), 1)+3*(st-1),data_comp{st}{sy,co}(:,1),'.','color',col_rep(1), 'markersize', 10) % Scatter plot of all values (dmPFC)
            hold on
            plot(repmat(2, size(data_comp{1}{sy,co},1), 1)+3*(st-1),data_comp{st}{sy,co}(:,2),'.','color',col_rep(2), 'markersize', 10) % Scatter plot of all values (dmPFC)
            plot(st+3*(st-1), nanmean(data_comp{st}{sy,co}), 'xk') % Plot mean
%             plot(repmat(1:2, size(data_comp{st}{sy,co},1),1)'+3*(st-1),data_comp{st}{sy,co}', 'color',col_rep(st)) % Plot individual rat traces
            set(gca, 'xlim', [0 6])
            xlabel('<-- s3 || s9 -->')

            % Test
            [~,pval(st),~,s_out]=ttest2(data_comp{st}{sy,co}(~isnan(data_comp{st}{sy,co}(:,1)),1),...
                data_comp{st}{sy,co}(~isnan(data_comp{st}{sy,co}(:,2)),2));
            test_outputs{st,sy}(:,co) = [round(s_out.tstat,4), s_out.df, round(pval(st),4)];
            if st == 2
                yl = get(gca, 'ylim');
                for gr2 = 1:2
                    if pval(gr2)<0.05
                        line([1 2] +3*(gr2-1), [ yl(2)+(yl(2)-yl(1))*0.1  yl(2)+(yl(2)-yl(1))*0.1], 'color', 'k')
                        text(1.5+3*(gr2-1),  yl(2)+(yl(2)-yl(1))*0.1, '*','HorizontalAlignment', 'center')
                        if pval(gr2)<0.01
                            text(1.5+3*(gr2-1),  yl(2)+(yl(2)-yl(1))*0.1, '**','HorizontalAlignment', 'center')
                            if pval(gr2)<0.001
                                text(1.5+3*(gr2-1),  yl(2)+(yl(2)-yl(1))*0.1, '***','HorizontalAlignment', 'center')
                            end
                        end
                    end
                end
                set(gca, 'ylim', [yl(1), yl(2)+(yl(2)-yl(1))*0.2])
            end
        end
        clear pval
    T_OUT_2{st,sy}=array2table(test_outputs{st,sy}, 'VariableNames', comps);
    end
    supertitle(snames{sy})
end


% 3.11e) Signal parameters per stage, comparisons
% Compares signal parameters across stages
nf=nf+2; %Figure indices
data_par = cell(1,1);
startst = ones(numel(comps),1)*3; % Earliest stage parameter is relevant (i.e. there is no ITI in stage 1 and 2, so dF during ITI isn't applicable)
endst = 10;

for sy = 1:2 % Loop syncs
    figure(nf+sy-1)
    for gr = 1:2 % Loop groups
        for co = 1:numel(comps)
            for st = startst(co):endst % Loop stages
                k=1;
                for ra = unique(meta_in(meta_in(:,1)==gr,2))' % Loop rats
                    trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st & meta_in(:,5)==1); % Select trials
                    if st == 10
                        trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st & meta_in(:,5)==1 & meta_in(:,6)==1); % Select trials
                    end
                    try % By-pass empty comparisons (sometimes it depends on sync whether or not there is one)
                        data_par{gr}{sy,co}(k,st)=nanmean(calcium(sy).(comps{co})(trid,:)); % Store mean values (Rat(X) * Stage(Y))
                    catch
                        data_par{gr}{sy,co}(k,st)=nan;
                    end
                    k=k+1;
                end
            end
            subplot(ceil(sqrt(numel(comps))), round(sqrt(numel(comps))), co) % Subplot
            
            % Remove NaN from dataset
            % Calcium values
            y1 = data_par{gr}{sy,co}(:,startst(co):end);
%             y2 = y1(~isnan(y1(:,1)) & ~isnan(y1(:,end)),:);
            y3 = y1(~isnan(y1(:)));
            
            % Corresponding x-axis values
            x1 = repmat(startst(co):endst, size(data_par{gr}{sy,co}(:,startst(co):end),1), 1)+10*(gr-1);
%             x2 = x1(~isnan(y1(:,1)) & ~isnan(y1(:,end)),:);
            x3 = x1(~isnan(y1(:)));
            
%             title([comps{co} ' - stage ' num2str(startst(co)) ' vs stage 9'])
%             plot(repmat(startst(co):endst, size(data_par{gr}{sy,co},1), 1)+endst+1*(gr-1),data_par{gr}{sy,co}(:,startst(co):end),'.','color',col_rep(gr), 'markersize', 10) % Scatter plot of all values
%             hold on
            
            % Generate trendline
%             coef = polyfit(x2,y2, 1);
            try
                lm = fitlm(x3(:),y3(:),'linear');
%                 % Generate 95% CIs
%                 p11 = predint(c2,x3,0.95,'observation','off');
%                 p12 = predint(c2,x3,0.95,'observation','on');
%                 p21 = predint(c2,x3,0.95,'functional','off');
%                 p22 = predint(c2,x3,0.95,'functional','on');
%                 
%                 plot(c2,x3,y3)%, 'color', col_rep(gr))
                plotAdded(lm)
                hold on
%                 plot(x3,p22,'m--') % Plot trendline and CIs
                legend('hide')
                title(comps{co})

                
                yl=get(gca,'ylim');
                text(median(startst(co):endst)+endst*(gr-1),yl(2)-0.05*diff(yl),['R2 = ' num2str(lm.Rsquared.Adjusted)]);
                text(median(startst(co):endst)+endst*(gr-1),yl(2)-0.15*diff(yl),['p = ' num2str(lm.Coefficients{2,4})]);
            catch 
                continue
            end
            set(gca, 'xlim',[startst(co)-1, 22], 'box', 'off')
%             plot(polyval(coef,startst(co):9)) % Plot trendline
            xlabel('<-- dmPFC || vmPFC -->')
        end
    end
    supertitle(snames{sy})
end


% 3.12) Means of rats, per stage, one figure per group.
% No distinction between outcome, trained status 
nf=nf+2; %Figure indices
cmap_st = ([.7:-.4/6:.3; 0:0.8/6:.8; 0.4:0.1/6:0.5])'; % Generate color map for stages. Purple = low, green = high

for gr = 1:2
    figure(nf+gr-1)
    for st = 3:9
        for sy = 1:2
            k=1;
            for ra = unique(meta_in(meta_in(:,1)==gr,2))'
                trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st & meta_in(:,5)==1);
                data_tmp{sy}(k,:)=nanmean(calcium(sy).traces(trid,:));
                k=k+1;
            end
            plot(xv{sy},movmean(nanmean(data_tmp{sy}),15), 'color', cmap_st(st-2,:), 'linewidth',3);
            hold on
            if st == 9 && sy==2
                yl = get(gca, 'ylim'); % Get y-axis limits
                line([5 5], yl, 'color','k')
                line([10 10], yl, 'color','k')
                line([21.5 21.5], yl, 'color','k')
                set(gca,'ylim',yl)
                set(gca,'box','off')
            end
        end
    end
    supertitle(['Means per stage - ' gnames{gr}])
    labels = findobj(gcf); %Find all line objects in figure
    legend(labels(end:-2:end-16), stnames) % Assign legend (just so happens to be that these specific line objects are the means. Not rly sure how this works exactly)
end

% 3.12a) Means of rats, stage by stage comparison
nf=nf+2; %Figure indices
plotcoord = cell(1,1);
plotcoord{1} = '00';
data_corr=cell(1,2);
dwin = {ceil(5*15.89):ceil(10*15.89), 1:ceil(8*15.89)}; % Correlation windows

for gr = 1:2
    n = 1;
    figure(nf+gr-1)
    for st = 3:9
        for st2 = 3:9
            if ismember([num2str(st) num2str(st2)],plotcoord) || st==st2 || ismember([num2str(st2) num2str(st)],plotcoord) % Skip if plot has already been made
                continue
            else
                subplot(9,9,st+9*(st2-1))
                for sy = 1:2
                    data_tmp=cell(1,2);
                    k=1;
                    for ra = unique(meta_in(meta_in(:,1)==gr,2))'
                        trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st);
                        data_tmp{1}(k,:)=nanmean(calcium(sy).traces(trid,:)); % Select data traces from 1 stage
                        trid2 = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st2);
                        data_tmp{2}(k,:)=nanmean(calcium(sy).traces(trid2,:)); % Select data traces from other stage
                        
                        [corr_tmp]=corrcoef(data_tmp{1}(k,dwin{sy}), data_tmp{2}(k,dwin{sy})); % Correlation coefficient (only rats with both stages)
                        data_corr{gr,sy}(st,st2,k)=corr_tmp(2,1); % Store correlation coefficients
                        data_corr{gr,sy}(st2,st,k)=corr_tmp(2,1);
                        data_corr{gr,sy}(st,st,k)=1;
                        k=k+1; % add number for rat loop
                    end
                    plot(xv{sy},nanmean(data_tmp{1}), 'color', col_rep(1), 'linewidth',2); % Plot means
                    hold on
                    jbfill(xv{sy}, nanmean(data_tmp{1})+nanstd(data_tmp{1})/sqrt(size(data_tmp{1},1)),...
                        nanmean(data_tmp{1})-nanstd(data_tmp{1})/sqrt(size(data_tmp{1},1)),...
                        col_rep(1), col_rep(1)); % Plot SEM
                    hold on
                    
                    plot(xv{sy},nanmean(data_tmp{2}), 'color', col_rep(2), 'linewidth',2); % Plot means
                    hold on
                    jbfill(xv{sy}, nanmean(data_tmp{2})+nanstd(data_tmp{2})/sqrt(size(data_tmp{2},1)),...
                        nanmean(data_tmp{2})-nanstd(data_tmp{2})/sqrt(size(data_tmp{2},1)),...
                        col_rep(2), col_rep(2)); % Plot SEM
                    hold on
                    
                    if sy == 2
                        yl = get(gca, 'ylim'); % Get y-axis limits
                        line([5 5], yl, 'color','k')
                        line([10 10], yl, 'color','k')
                        line([21.5 21.5], yl, 'color','k')
                        set(gca,'ylim',yl)
                    end
                    
                    if st == 1
                        ylabel(['Stage ' num2str(st2)])
                    end
                    if st2 == 9
                        xlabel(['Stage ' num2str(st)])
                    end
                end
            end
            plotcoord{n}=[num2str(st) num2str(st2)];
            n=n+1;
       end
    end
    supertitle(gnames{gr})
    plotcoord = cell(1,1);
    plotcoord{1} = '00';
end

% 3.12b) Means of rats, correlations between stages
nf=nf+2; %Figure indices
splotpos = {10:13, 15:18}; %subplot positions

for gr = 1:2
    figure(nf+gr-1)
    for sy2 = 1:2 % Loop synchronizations
        for pp = 3:9 % Loop stages
            subplot(2,9,pp) %Subplot
            dcy = reshape(permute(data_corr{gr,sy2}(pp,:,:),[2,1,3]),size(data_corr{gr,sy2}(pp,:,:),2),[])'; % Retrieve correlation coefficients per rat, per stage
            dcx = repmat((1:9)+10*(sy2-1),size(dcy,1),1); % x-valules for plot
            dcmean = nanmean(dcy); %means
            dcsemup = dcmean+nanstd(dcy)/sqrt(size(dcy,1)); %SEM (upper)
            dcsemdown = dcmean-nanstd(dcy)/sqrt(size(dcy,1)); % SEM (lower)
            plot(dcx,dcy,'.', 'color', col_rep(1), 'markersize',6) %plot individual rats
            hold on
            plot((1:9)+10*(sy2-1), dcmean, 'x', 'color', 'k', 'markersize', 8) % plot means ('X')
            jbfill((1:9)+10*(sy2-1), dcsemup, dcsemdown,col_rep(1),col_rep(1)); % plot SEM
            hold on
            
            if pp==5
                title(['Correlation coefficients per stage' newline gnames{gr}])
            end
            xlabel(['Correlation to stage ' num2str(pp)])
        end
        
        subplot(2,9,splotpos{sy2}) %Subplot
        ddx3 = 0:8;
        for p = 2:9
            ddx3(p,:) =  ddx3(p-1,:)-1;
        end
        ddx = repmat(ddx3,1,1,size(data_corr{gr,sy2},3));
        
        for m = 1:8
            dz = data_corr{gr,sy2}(ddx == m);
            dzsemup(m) = nanmedian(dz)+nanstd(dz)/sqrt(numel(dz));
            dzsemdown(m) = nanmedian(dz)-nanstd(dz)/sqrt(numel(dz));
            rgrdz{gr,sy2}{m}=[dz, m*ones(numel(dz),1)];
            
            plot(m*ones(numel(dz),1), dz, '.', 'color', col_rep(1), 'markersize', 8)
            hold on
            plot(m,nanmedian(dz), 'xk', 'markersize', 10)
            
            if m == 8
                jbfill(1:8, dzsemup, dzsemdown, col_rep(1), col_rep(1),1, 0.35);
                hold on
                xlim([0 9])
            end
        end
        title(['Correlation to X stages apart' newline snames{sy2}])
    end
end

% 3.13) Means of rats, per stage, separate figures, outcomes split
nf=nf+2; %Figure indices
lw = [2 1 1 1]; %linewidth
pw = [0.4 0.2 0.2 0.2]; %SEM patch transparency
onames = {'Correct' 'Incorrect' 'Omission' 'Premature'};

for gr = 1:2
    figure(nf+gr-1)
    for st = 1:9
        subplot(3,3,st)
        for oc = 1:4
            data_tmp = cell(1);
            for sy = 1:2
                k=1;
                for ra = unique(meta_in(meta_in(:,1)==gr,2))'
                    trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st & meta_in(:,5)==oc);
                    if sum(trid)<=1
                        data_tmp{sy}(k,:)=nan(1,numel(xv{sy}));
                    else
                        data_tmp{sy}(k,:)=nanmean(calcium(sy).traces(trid,:));
                    end
                    k=k+1;
                end
                p(oc) = plot(xv{sy},nanmean(data_tmp{sy}), 'color', col_resp(oc), 'linewidth',lw(oc)); % Plot means
                hold on
                jbfill(xv{sy}, nanmean(data_tmp{sy})+nanstd(data_tmp{sy})/sqrt(size(data_tmp{sy},1)),...
                    nanmean(data_tmp{sy})-nanstd(data_tmp{sy})/sqrt(size(data_tmp{sy},1)),...
                    col_resp(oc), col_resp(oc),1,pw(oc)); % Plot SEM
                hold on
                title(['stage ' num2str(st)])
                if sy == 2 && oc == 4
                    yl = get(gca, 'ylim'); % Get y-axis limits
                    line([5 5], yl, 'color','k')
                    line([10 10], yl, 'color','k')
                    line([21.5 21.5], yl, 'color','k')
                    set(gca,'ylim',yl)
                    if st == 1
                        legend(p, onames, 'location', 'south', 'orientation', 'horizontal')
                    end
               end
            end
        end
    end
    supertitle(gnames{gr})
end

% 3.14) Means of rats, per stage, one figure per group, split by outcomes
nf=nf+2; %Figure indices
cmap_st = ([.7:-.4/8:.3; 0:0.8/8:.8; 0.4:0.1/8:0.5])'; % Generate color map for stages. Purple = low, green = high
% 
for gr = 1:2
    figure(nf+gr-1)
    for st = 1:9
        for oc = 1:4
            data_tmp = cell(1);
            subplot(2,2,oc)
            title(onames{oc})
            for sy = 1:2
                k=1;
                for ra = unique(meta_in(meta_in(:,1)==gr,2))'
                    trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st & meta_in(:,5)==oc);
                    if sum(trid)<=1
                        data_tmp{sy}(k,:)=nan(1,numel(xv{sy}));
                    else
                        data_tmp{sy}(k,:)=nanmean(calcium(sy).traces(trid,:));
                    end
                    k=k+1;
                end
                p(oc) = plot(xv{sy},nanmean(data_tmp{sy}), 'color', cmap_st(st,:), 'linewidth',2);
                hold on
                if st == 9
                    yl = get(gca, 'ylim'); % Get y-axis limits
                    line([5 5], yl, 'color','k')
                    line([10 10], yl, 'color','k')
                    line([21.5 21.5], yl, 'color','k')
                    set(gca,'ylim',yl)
                end
            end
        end
        supertitle(gnames{gr})
    end
end

% 3.15) Means of rats, per stage, separate figures, outcomes AND trained
% status split
% ATM only for correct responses and all errors pooled
nf=nf+2; %Figure indices
lw = [2 1 1 1]; % linewidth
pw = [0.4 0.2 0.2 0.2]; %transparency of SEM area
ls = {'-', ':'}; % linespec
ctype = [1 1; 2 2; 3 3; 4 4];

for gr = 1:2
    figure(nf+gr-1)
    for st = 1:9
        subplot(3,3,st)
        for oc = 1:4
            data_tmp = cell(1);
            for sy = 1:2
                k=1;
                for ra = unique(meta_in(meta_in(:,1)==gr,2))'
                    trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st & ismember(meta_in(:,5),ctype(oc,:)) & meta_in(:,11)==0);
                    trid2 = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st & ismember(meta_in(:,5),ctype(oc,:)) & meta_in(:,11)==1);
                    if sum(trid)<=1
                        data_tmp{sy,1}(k,:)=nan(1,numel(xv{sy}));
                        if sum(trid2)<=1
                            data_tmp{sy,2}(k,:)=nan(1,numel(xv{sy}));
                        end
                    elseif sum(trid2)<=1
                            data_tmp{sy,2}(k,:)=nan(1,numel(xv{sy}));
                    else
                        data_tmp{sy,1}(k,:)=nanmean(calcium(sy).traces(trid,:));
                        data_tmp{sy,2}(k,:)=nanmean(calcium(sy).traces(trid2,:));
                    end
                    k=k+1;
                end
                % UNtrained
                putr(oc) = plot(xv{sy},movmean(nanmean(data_tmp{sy,1}),15), 'color', col_resp(oc), 'linewidth',lw(oc), 'linestyle', ls{1}); % Plot means
                jbfill(xv{sy}, nanmean(data_tmp{sy})+nanstd(data_tmp{sy})/sqrt(size(data_tmp{sy},1)),...
                    nanmean(data_tmp{sy})-nanstd(data_tmp{sy})/sqrt(size(data_tmp{sy},1)),...
                    col_resp(oc), col_resp(oc),1,pw(oc)); % Plot SEM 
                hold on
                
                % Trained
                ptr(oc) = plot(xv{sy},movmean(nanmean(data_tmp{sy,2}),15), 'color', col_resp(oc), 'linewidth',lw(oc), 'linestyle', ls{2}); % Plot means
                hold on
                title(['stage ' num2str(st)])
                if sy == 2 && oc == 2
                    yl = get(gca, 'ylim'); % Get y-axis limits
                    line([5 5], yl, 'color','k')
                    line([10 10], yl, 'color','k')
                    line([21.5 21.5], yl, 'color','k')
                    set(gca,'ylim',yl)
%                     if st == 1
%                         legend([putr ptr], {'Correct (untr.)', 'Error (untr.)', 'Correct (tr.)', 'Error (tr.)'}, 'location', 'south', 'orientation', 'horizontal')
%                     end
               end
            end
        end
    end
    supertitle([gnames{gr} newline])
end


% 3.15b) Means of rats, per stage, separate figures, outcomes AND trained
% status split
% With permutation test, only stage 9
% This is rly mostly to find differences between error and corrects
nf=nf+2; %Figure indices
lw = [2 1 1 1]; % linewidth
% pw = [0.4 0.2 0.2 0.2]; %transparency of SEM area
ls = {'-', ':'}; % linespec
% ctype = [1 1; 2 2; 3 3; 4 4];
comps = {'auc' 'maxA' 'maxT' 'amplitude' 'dfiti' 'fcue' 'meanpeakwidth' 'meanpeakprominence' 'meanpeakauc' ...
    'weight' 'weight1' 'weight2' 'halfweighttime' 'fresp' 'peakrespdf' 'respdecaydf', 'variance'}; % Calcium signal parameters
% cst=3; %current stage

for gr = 1:2
    data_tmp_comp = cell(1,1);
    for sy = 1:2
        figure(nf+sy-1)
        for co = 1:numel(comps)
            for oc = 1:4
                k=1;
                for ra = unique(meta_in(meta_in(:,1)==gr,2))'% Loop rats
                    trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==3 & meta_in(:,5)==oc & meta_in(:,11)==0);
                    trid2 = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==3 & meta_in(:,5)==oc & meta_in(:,11)==1);
                    
                    try % By-pass empty comparisons (sometimes it depends on sync whether or not there is one)
                        data_tmp_comp{sy,co}(1,k,oc)= ...
                            nanmean(calcium(sy).(comps{co})(trid,:)); % Store mean values for each stage
                    catch
                        data_tmp_comp{sy,co}(1,k,oc)=nan(1,1);
                    end
                    
                    try % By-pass empty comparisons (sometimes it depends on sync whether or not there is one)
                        data_tmp_comp{sy,co}(2,k,oc)= ...
                            nanmean(calcium(sy).(comps{co})(trid2,:)); % Store mean values for each stage
                    catch
                        data_tmp_comp{sy,co}(2,k,oc)=nan(1,1);
                    end
                    k=k+1;
                end
                data_tmp_comp{sy,co}(data_tmp_comp{sy,co}==0)=nan;
                
                
                subplot(ceil(sqrt(numel(comps))), floor(sqrt(numel(comps))), co) % Subplot
                title(comps{co})
                plot(repmat((1:2)+2*(oc-1), size(data_tmp_comp{sy,co},2), 1)',data_tmp_comp{sy,co}(:,:,oc), 'color','k') % Plot individual rat traces
                hold on
                plot(repmat((1:2)+2*(oc-1), size(data_tmp_comp{sy,co},2), 1),data_tmp_comp{sy,co}(:,:,oc)','.','color',col_resp(oc), 'markersize', 12) % Scatter plot of all values
                %                 plot((1:2), nanmean(data_tmp_comp{sy,co}), 'xk') % Plot mean
                set(gca, 'xlim', [0 9])
                title(comps{co})
            end
        end
    supertitle([gnames{gr} ' ' snames{sy} newline])
    end

end

% 3.16a) Means of rats, per stage, separate figures, dmPFC vs vmPFC
tnames = {'Untrained' 'All'};
nf=nf+2; %Figure indices

for trn = 1:2
    figure(nf+trn-1)
    for st = 3
        data_tmp = cell(1);
        for gr = 1:2
%             subplot(3,3,st)
            for oc = 1 % trial OutCome
                for sy = 1:2
                    k=1;
                    for ra = unique(meta_in(meta_in(:,1)==gr,2))'
                        trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st & meta_in(:,5)==oc & meta_in(:,11)<trn);
                        if sum(trid)<=1
                            data_tmp{gr,sy}(k,:)=nan(1,numel(xv{sy}));
                        else
                            data_tmp{gr,sy}(k,:)=nanmean(calcium(sy).traces(trid,:));
                        end
                        k=k+1;
                    end
                    p2(gr) = plot(xv{sy},movmean(nanmean(data_tmp{gr,sy}),15), 'color', col_rep(gr), 'linewidth',lw(oc), 'linestyle', ls{1}); % Plot means
                    jbfill(xv{sy}, movmean(nanmean(data_tmp{gr,sy})+nanstd(data_tmp{gr,sy})/sqrt(size(data_tmp{gr,sy},1)),15),...
                        movmean(nanmean(data_tmp{gr,sy})-nanstd(data_tmp{gr,sy})/sqrt(size(data_tmp{gr,sy},1)),15),...
                        col_rep(gr), col_rep(gr)); % Plot SEM
                    hold on
                    
                    yl = get(gca, 'ylim'); % Get y-axis limits
                    
                    % Permutation test
                    if gr == 2
                        sig_out = permTest_array(data_tmp{1,sy}(~isnan(data_tmp{1,sy}(:,1)),:),...
                            data_tmp{2,sy}(~isnan(data_tmp{2,sy}(:,1)),:),5000);
                        xw = xv{sy}; xw(sig_out>0.05)=nan;
                        plot(xw, ones(numel(xw),1)*yl(2)*0.95, 'linewidth', 3, 'color', 'k')


%                         plot(xv{sy}(sig_out<0.05), ones(sum(sig_out<0.05),1)*yl(2)*0.95, 'linewidth', 3, 'color', 'k')
                        
                    end
                    
                    title(['stage ' num2str(st)])
                    
                    if sy == 2
                        line([5 5], yl, 'color','k')
                        line([10 10], yl, 'color','k')
                        line([21.5 21.5], yl, 'color','k')
                        set(gca,'ylim',yl)
                        if st == 3 && gr == 2
                            legend(p2(1:2), gnames, 'location', 'south', 'orientation', 'horizontal')
                        end
                        set(gca, 'box','off')
                        ylabel('dF/F')
                    end
                end
            end
        end
    end
    supertitle(['vmPFC and dmPFC comparison - ' tnames{trn} newline])
end

% 3.16b) Signal parameters, dmPFC vs vmPFC
nf=nf+2; %Figure indices

comps = {'auc' 'maxA' 'maxT' 'amplitude' 'dfiti' 'fcue' 'meanpeakwidth' 'meanpeakprominence' 'meanpeakauc' ...
    'weight' 'weight1' 'weight2' 'halfweighttime' 'fresp' 'peakrespdf' 'respdecaydf', 'variance'}; % Calcium signal parameters
data_comp_gr = cell(1,1);
cst=3; %current stage

for sy = 1:2 % Loop syncs
    figure(nf+sy-1)
    for co = 1:numel(comps)
        pval = ones(1,2);
        % Compare early and late stage trial data for each parameter
        for gr = 1:2
            k=1;
            for ra = unique(meta_in(meta_in(:,1)==gr,2))'% Loop rats
                trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==cst & meta_in(:,11)<=1); % Select trials
                
                try % By-pass empty comparisons (sometimes it depends on sync whether or not there is one)
                    data_comp_gr{sy,co}(k,gr)= ...
                        nanmean(calcium(sy).(comps{co})(trid,:)); % Store mean values for each stage
                catch
                    data_comp_gr{sy,co}(k,gr)=nan(1,1);
                end
                k=k+1;
            end
        end
        
        data_comp_gr{sy,co}(data_comp_gr{sy,co}==0)=nan;
        
        
        subplot(ceil(sqrt(numel(comps))), floor(sqrt(numel(comps))), co) % Subplot
        title(comps{co})
        plot(repmat(1:2, size(data_comp_gr{sy,co},1), 1),data_comp_gr{sy,co},'.','color','k', 'markersize', 12) % Scatter plot of all values
        hold on
        plot((1:2), nanmean(data_comp_gr{sy,co}), 'xk') % Plot mean
%         plot(repmat(1:2, size(data_comp_gr{sy,co},1),1)',data_comp_gr{sy,co}', 'color','k') % Plot individual rat traces
        set(gca, 'xlim', [0 3])
        title(comps{co})
        
        % Test
        [~,pval,~,~]=ttest2(data_comp_gr{sy,co}(:,1), data_comp_gr{sy,co}(:,2));
        yl = get(gca, 'ylim');
        if pval<0.05
            line([1 2], [ yl(2)+(yl(2)-yl(1))*0.1  yl(2)+(yl(2)-yl(1))*0.1], 'color', 'k')
            text(1.5,  yl(2)+(yl(2)-yl(1))*0.1, '*','HorizontalAlignment', 'center')
            if pval<0.01
                text(1.5,  yl(2)+(yl(2)-yl(1))*0.1, '**','HorizontalAlignment', 'center')
                if pval<0.001
                    text(1.5,  yl(2)+(yl(2)-yl(1))*0.1, '***','HorizontalAlignment', 'center')
                end
            end
        end
        set(gca, 'ylim', [yl(1), yl(2)+(yl(2)-yl(1))*0.2])
        
        clear s_out
        supertitle([snames{sy} ' - stage ' num2str(cst) newline])
    end
end

% 3.17a) How does signal develop across stages? Similar figure to behavior
% figure 2.6a
% 6. How does performance develop across stages?
% a. Performance parameters across all trials
nf=nf+1; %Figure indices
comptype = 'accwin'; % Define behavioral parameter
sigtype = 'fcue'; % Define signal parameter (type 'calcium' in command window for full list')
dz = cell(1,1); ds = cell(1,1); % Initiate array
lspec = {'-', ':'}; % Line style in plot
bincount = 20; % Number of bins in sliding window

figure(nf)
for gr = 1:2
    for st = 3:9 %unique(behavior.singleTrials{gr}{ra,trn}(:,3))'
        for ra = 1:size(behavior.accwin{gr},1)
            subplot(3,2,gr)
            
            % Plot accuracy for each individual rat
            stcount = sum(behavior.singleTrials{gr}{ra,1}(:,3)==st); % Trial count in stage
            plot(st:1/(stcount-1):st+1, behavior.(comptype){gr}{ra,1}(behavior.singleTrials{gr}{ra,1}(:,3)==st)+ra*1.33,...
                'color', col_rep(st),'linewidth',2, 'linestyle', lspec{1})
            hold on
            title(['Individual rats ' comptype ' development in ' gnames{gr}])
            ylabel('Rats')
            xlabel('Stages')
            
            % Plot sliding window of accuracy for each rat, scaled to
            % trained and untrained phases of task
            dt = st:1/(stcount-1):st+1; % x-axis for trial count
            dy = behavior.(comptype){gr}{ra,1}(behavior.singleTrials{gr}{ra,1}(:,3)==st); % performance values 
            [~,~,binid]=histcounts(dt,bincount); % Assign trials to bins (if more than 20)

            % Calculate sliding window of signal parameters in similar
            % fashion
            rat_id = unique(meta(meta(:,1)==gr,2));
            trid = (meta(:,1)==gr & meta(:,2)==rat_id(ra) & meta(:,3)==st & meta(:,5)==1);
            dsig = calcium(1).(sigtype)(trid,:);
            dsigt = st:1/(sum(trid)-1):st+1; % x-axis for trial count
            [~,~,binid2]=histcounts(dsigt,bincount); % Assign trials to bins (if more than 20)
%             dsig = nansum(dsig,2);
            
            % Assign values to bins
            if isempty(binid) % If no trials, assign NaN value
                dz{gr}{st,1}(ra,:) = nan(1,bincount);
            else
                for z = 1:max(binid) % Loop through bins to assign mean performance value
                    dz{gr}{st,1}(ra,z) =  mean(dy(binid==z));
                end
            end
            if isempty(binid2) % If no trials, assign NaN value
                ds{gr}{st,1}(ra,:) = nan(1,bincount);
            else
                for z = 1:max(binid2) % Loop through bins to assign mean performance value
                    ds{gr}{st,1}(ra,z) =  mean(dsig(binid2==z));
                end
            end
        end

        % Plot behavioral timeline
        subplot(3,2,gr+2) 
        dz{gr}{st,1}(dz{gr}{st,1}==0)=nan; % NaN if no trials
        plot(st:1/(bincount-1):st+1, nanmean(dz{gr}{st,1}),...
            'color', col_rep(st),'linewidth',2, 'linestyle', lspec{1}) % Plot bin means
        hold on
        jbfill(st:1/(bincount-1):st+1,...
            nanmean(dz{gr}{st,1})+nanstd(dz{gr}{st,1})/sqrt(size(dz{gr}{st,1},1)),...
            nanmean(dz{gr}{st,1})-nanstd(dz{gr}{st,1})/sqrt(size(dz{gr}{st,1},1)),...
            col_rep(st), col_rep(st)); % Plot SEM
        hold on
        title(['Development of ' comptype ' in ' gnames{gr}])
        ylabel('%')
        xlabel('Stages')
        
        
        % Plot calcium signal parameter timeline
        subplot(3,2,gr+4)
        ds{gr}{st,1}(ds{gr}{st,1}==0)=nan; % NaN if no trials
        plot(st:1/(bincount-1):st+1, nanmean(ds{gr}{st,1}),...
            'color', col_rep(st),'linewidth',2, 'linestyle', lspec{1}) % Plot bin means
        hold on
        jbfill(st:1/(bincount-1):st+1,...
            nanmean(ds{gr}{st,1})+nanstd(ds{gr}{st,1})/sqrt(size(ds{gr}{st,1},1)),...
            nanmean(ds{gr}{st,1})-nanstd(ds{gr}{st,1})/sqrt(size(ds{gr}{st,1},1)),...
            col_rep(st), col_rep(st)); % Plot SEM
        hold on
        title(['Development of ' sigtype ' in ' gnames{gr}])
        xlabel('Stages')
    end
end


% 3.17b). Distinction between trained and untrained trials. Similar to
% 3.17a, but now trained and untrained trials are separated
nf=nf+1; %Figure indices
dz = cell(1,1); ds = cell(1,1);

figure(nf)
for gr = 1:2 % Loop through groups (1=only dmPFC, 2=only vmPFC, 3=all)
    for trn=1:2 % Loop through training status (untrained / trained)
        for st = 1:9 %unique(behavior.singleTrials{gr}{ra,trn}(:,3))' % Loop through stages (1= stage 1, 11= vITI)
            for ra = 1:size(behavior.accwin{gr},1) % Loop through rats
                subplot(3,2,gr)
                stcount = sum(behavior.singleTrials{gr}{ra,trn+1}(:,3)==st); % trialcount
                plot(st+.5*(trn-1):.5/(stcount-1):st+.5*trn, behavior.(comptype){gr}{ra,trn+1}(behavior.singleTrials{gr}{ra,trn+1}(:,3)==st)+ra*1.33,...
                    'color', col_rep(st),'linewidth',2, 'linestyle', lspec{trn}) % Plot individual rat sliding performance window
                hold on
                title(['Individual rats ' comptype ' development in ' gnames{gr}])
                ylabel('Rats')
                xlabel('Stages')
                
                % categorize into bins
                dt = st:1*trn/(stcount-1):st+1; % x-axis values (so that trained and untrained sections each comprise 0.5 units on x-axis
                dy = behavior.(comptype){gr}{ra,trn+1}(behavior.singleTrials{gr}{ra,trn+1}(:,3)==st); % Behavioral parameter values
                [~,~,binid]=histcounts(dt,bincount); % Divide trials into bins
                
                % Calculate sliding window of signal parameters in similar
                % fashion
                rat_id = unique(meta(meta(:,1)==gr,2));
                trid = (meta(:,1)==gr & meta(:,2)==rat_id(ra) & meta(:,3)==st & meta(:,5)==1 & meta(:,11)==trn-1);
                dsig = calcium(1).(sigtype)(trid,:);
                %                 dsig = nansum(dsig,2);
                dsigt = st:1/(sum(trid)-1):st+1; % x-axis for trial count
                [~,~,binid2]=histcounts(dsigt,bincount); % Assign trials to bins (if more than 20)
                
                % Assign values to bins
                if isempty(binid) % If no trials, assign NaN value
                    dz{gr}{st,trn}(ra,:) = nan(1,bincount);
                else
                    for z = 1:max(binid) % Loop through bins to assign mean performance value
                        dz{gr}{st,trn}(ra,z) =  mean(dy(binid==z));
                    end
                end
                if isempty(binid2) % If no trials, assign NaN value
                    ds{gr}{st,trn}(ra,:) = nan(1,bincount);
                else
                    for z = 1:max(binid2) % Loop through bins to assign mean performance value
                        ds{gr}{st,trn}(ra,z) =  mean(dsig(binid2==z));
                    end
                end
            end
            
            % Plot behavioral timeline
            subplot(3,2,gr+2)
            dz{gr}{st,trn}(dz{gr}{st,trn}==0)=nan;
            plot(st+.5*(trn-1):.5/(bincount-1):st+.5*trn, nanmean(dz{gr}{st,trn}),...
                'color', col_rep(st),'linewidth',2, 'linestyle', lspec{trn}) % Plot mean traces
            hold on
            jbfill(st+.5*(trn-1):.5/(bincount-1):st+.5*trn,...
                nanmean(dz{gr}{st,trn})+nanstd(dz{gr}{st,trn})/sqrt(size(dz{gr}{st,trn},1)),...
                nanmean(dz{gr}{st,trn})-nanstd(dz{gr}{st,trn})/sqrt(size(dz{gr}{st,trn},1)),...
                col_rep(st), col_rep(st)); % SEM confidence intervals
            hold on
            xlim([0 10])
            
            % Plot calcium signal parameter timeline
            subplot(3,2,gr+4)
            ds{gr}{st,trn}(ds{gr}{st,trn}==0)=nan; % NaN if no trials
            plot(st+.5*(trn-1):.5/(bincount-1):st+.5*trn, nanmean(ds{gr}{st,trn}),...
                'color', col_rep(st),'linewidth',2, 'linestyle', lspec{trn}) % Plot mean traces
            hold on
            jbfill(st+.5*(trn-1):.5/(bincount-1):st+.5*trn,...
                nanmean(ds{gr}{st,trn})+nanstd(ds{gr}{st,trn})/sqrt(size(ds{gr}{st,trn},1)),...
                nanmean(ds{gr}{st,trn})-nanstd(ds{gr}{st,trn})/sqrt(size(ds{gr}{st,trn},1)),...
                col_rep(st), col_rep(st)); % SEM confidence intervals
            hold on
            title(['Development of ' sigtype ' in ' gnames{gr}])
            xlabel('Stages')
             xlim([0 10])
        end
    end
end


% Save figures
figHandles = findobj('Type', 'figure'); % Gather all figure handles
set(figHandles, 'renderer', 'painters'); % Set renderer to 'painters' for proper further processing in Illustrator
set(figHandles, 'Position', get(0, 'Screensize')); % Set figures to full screen size (easier to process images in Illustrator)

for fign = 1:numel(figHandles)
    saveas(figHandles(fign), [savefolder2, '\Calcium-',num2str(figHandles(fign).Number)], 'epsc'); % Save each figure as .EPSC file 
end


%% Variable ITI and SD
% Meta table contents (useful when tweaking this chapter)
% Generate or append trial data matrix
% 1      2    3     4     5     6        7           8
% group, rat, cond, sess, resp, subspec, lower norm, upper norm
% 9             % 10         % 11
% trial time,   resp time    learning/nolearning

% 3.20a) Compare signal between ITI and SD conditions within group
xv{1} = 0:17.5/(size(calcium(1).traces,2)-1):17.5; % x-axis for combined sync plots (trialstart
xv{2} = 18.5+(0:8/(size(calcium(2).traces,2)-1):8);  % x-axis for combined sync plots (response)
nf = 1; % Figure number
ls = {'-' '-.' ':'};
ocspec = {1, 2};
cname320 = {'Correct' 'Errors'};
sname320 = {'vSD' 'vITI'};
gnames = {'dmPFC' 'vmPFC'};
% scname320a = {'variable ITI/SD', 'fixed ITI/SD'};

for gr = 1:2
    figure(nf+gr-1)
    traces_vs = cell(1,1);
    for st = 10:11 % Loop through vITI and vSD
        for oc = 1:2 % OutCome (corr/inc etc)
            for sy = 1:2
                
                % Data selection
                data_in = data_in_z{sy}; % Non-proportionalized
                subplot(2,2,2*(st-10)+oc)
                k = 1;
                for ra = unique(meta_in(meta_in(:,1)==gr,2))'
                    trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st & meta_in(:,6)>=1 & meta_in(:,5)==ocspec{oc}); % viti/sd
                    trid2 = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==9 & meta_in(:,5)==ocspec{oc}); %s9
                    
                    
                    if sum(sum(trid,2))>5
                        traces_vs{st-8,sy}{oc}(k,:) = nanmean(data_in(sum(trid,2)==1,:));
                    else
                        traces_vs{st-8,sy}{oc}(k,:) = nan(1,size(data_in,2));
                    end
                    
                    if sum(sum(trid2,2))>5
                        traces_vs{1,sy}{oc}(k,:) = nanmean(data_in(sum(trid2,2)==1,:)); % traces_vs are rat mean traces in either stage 9 or vITI/vSD
                    else
                        traces_vs{1,sy}{oc}(k,:) = nan(1,size(data_in,2));
                    end

                    k = k+1;
                end
                
                
                % Plots
                s(1)=plot(xv{sy},nanmean(traces_vs{st-8,sy}{oc}), 'color', col_resp(oc), 'linewidth', 2); %Plot vSD/vITI means
                hold on
                jbfill(xv{sy}, nanmean(traces_vs{st-8,sy}{oc})+nanstd(traces_vs{st-8,sy}{oc})/sqrt(size(traces_vs{st-8,sy}{oc},1)),...
                    nanmean(traces_vs{st-8,sy}{oc})-nanstd(traces_vs{st-8,sy}{oc})/sqrt(size(traces_vs{st-8,sy}{oc},1)),...
                    col_resp(oc), col_resp(oc)); % Plot SEM
                hold on
                
                s(2)=plot(xv{sy},nanmean(traces_vs{1,sy}{oc}), 'color', 'k', 'linewidth', 1); %Plot S9 means
                jbfill(xv{sy}, nanmean(traces_vs{1,sy}{oc})+nanstd(traces_vs{1,sy}{oc})/sqrt(size(traces_vs{1,sy}{oc},1)),...
                    nanmean(traces_vs{1,sy}{oc})-nanstd(traces_vs{1,sy}{oc})/sqrt(size(traces_vs{1,sy}{oc},1)),...
                    'k', 'k'); % Plot SEM
                hold on
                
                % Permutation test
                sig_out = permTest_array(traces_vs{st-8,sy}{oc}(~isnan(traces_vs{st-8,sy}{oc}(:,1)),:),...
                    traces_vs{1,sy}{oc}(~isnan(traces_vs{1,sy}{oc}(:,1)),:),5000);
                xw = xv{sy}; xw(sig_out>0.05)=nan;
                yl = get(gca, 'ylim');
                plot(xw, ones(numel(xw),1)*yl(2)*0.95, 'linewidth', 3, 'color', 'k')

                % Plot parameters
                if sy == 2
                    title([sname320{st-9} ' - ' cname320{oc}])
                    ylabel('dF/F')
                    xlabel('Time')
                    line([5 5], yl, 'color', 'k')
                    line([10 10], yl, 'color', 'k')
                    line([21.5 21.5], yl, 'color', 'k')
                    if st == 11
                        line([12.5 12.5], yl, 'color','k', 'linestyle', ls{2})
                        line([17.5 17.5], yl, 'color', 'k', 'linestyle',ls{3})
                    end
                    set(gca, 'ylim', yl)
                    legend(s, {sname320{st-9} 'fixedITI'}, 'location', 'best', 'orientation', 'vertical')
                end
            end
        end
    end
    supertitle([gnames{gr} newline])
end


% 3.20b) Compare signal between ITI and SD conditions within group,
% subconditions separated
nf = nf+2; % Figure number
ls = {'-' '-.' ':'};
ocspec = {1, [2 3 4]};
cname320 = {'Correct' 'Errors'};
sname320 = {'vSD' 'vITI'};
scname320b = {'easy', 'mid', 'hard', 'fixed ITI'};

for gr = 1:2
    figure(nf+gr-1)
    traces_vs = cell(1,1);
    for sy = 1:2
        data_subspecs = data_in_z{sy}; % Non-proportionalized
        for st = 10:11
            for oc = 1:2 % OutCome (corr/inc etc)
                subplot(2,2,2*(st-10)+oc)
                for sc = 1:3 % SubCondition (ITIs or SDs)
                    k = 1;
                    for ra = unique(meta_in(meta_in(:,1)==gr,2))'
                        trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st & meta_in(:,5)==ocspec{oc} & meta_in(:,6)==sc); % variable sd/iti
                        trid2 = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==9 & meta_in(:,5)==ocspec{oc}); % stage 9
                        
                        
                        if sum(sum(trid,2))> 5 % Check if more than 5 trials in sessions for each rat
                            traces_vs{st-8,sy}{oc}(k,:) = nanmean(data_subspecs(sum(trid,2)==1,:)); % vsd/iti
                        else
                            traces_vs{st-8,sy}{oc}(k,:) = nan(1,size(data_subspecs,2)); %vsd/iti
                        end
                        if sum(sum(trid2,2))>5 % Check if more than 5 trials in sessions for each rat
                            traces_vs{1,sy}{oc}(k,:) = nanmean(data_subspecs(sum(trid2,2)==1,:)); % s9
                        else
                            traces_vs{1,sy}{oc}(k,:) = nan(1,size(data_subspecs,2)); % s9
                        end
                        k = k+1;
                    end
                    
                    r(sc)=plot(xv{sy},nanmean(traces_vs{st-8,sy}{oc}), 'color', col_resp(oc), 'linestyle', ls{sc}, 'linewidth', 2);
                    hold on
                    r(4)=plot(xv{sy},nanmean(traces_vs{1,sy}{oc}), 'color', 'k', 'linewidth', 1);
                end
                
                title([sname320{st-9} ' - ' cname320{oc}])
                ylabel('dF/F')
                xlabel('Time')
                yl = get(gca, 'ylim');
                line([5 5], yl, 'color', 'k')
                line([10 10], yl, 'color', 'k')
                line([22.5 22.5], yl, 'color', 'k')
                if st == 11
                    line([12.5 12.5], yl, 'color','k', 'linestyle', ls{2})
                    line([17.5 17.5], yl, 'color', 'k', 'linestyle',ls{3})
                end
                set(gca, 'ylim', yl)
                legend(r, scname320b, 'location', 'best', 'orientation', 'vertical')
            end
        end
    end
    supertitle([gnames{gr} newline])
end


% 3.21a) Compare signal between vITI and vSD conditions within group
xv{1} = 0:17.5/(size(calcium(1).traces,2)-1):17.5; % x-axis for combined sync plots (trialstart
xv{2} = 18.5+(0:8/(size(calcium(2).traces,2)-1):8);  % x-axis for combined sync plots (response)
nf = nf+1; % Figure number
lspec = {'-' '-.' ':'};
ocspec = {1, 2};
cname320 = {'Correct' 'Errors'};
sname320 = {'vSD' 'vITI'};
gnames = {'dmPFC' 'vmPFC'};
% scname320a = {'variable ITI/SD', 'fixed ITI/SD'};
stnos = [9 10];
compo = [1 2];

figure(nf)

for gr = 1:2
    traces_vs = cell(1,1);
    for st = 1:numel(stnos) % Loop through stage types
        for sy = 1:2
            
            % Data selection
            data_in = data_in_z{sy}; % Non-proportionalized
            subplot(1,2,gr)
            k = 1;
            for ra = unique(meta_in(meta_in(:,1)==gr,2))'
                trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==stnos(st) & meta_in(:,6)<=1 & meta_in(:,5)==1); % viti/sd
                
                if sum(sum(trid,2))>5
                    traces_vs{st,sy}(k,:) = nanmean(data_in(sum(trid,2)==1,:));
                else
                    traces_vs{st,sy}(k,:) = nan(1,size(data_in,2));
                end
                k = k+1;
            end
            
            
            % Plots
            s(st)=plot(xv{sy},nanmean(traces_vs{st,sy}), 'color', col_resp(st), 'linewidth', 2); %Plot vSD/vITI means
            hold on
            jbfill(xv{sy}, nanmean(traces_vs{st,sy})+nanstd(traces_vs{st,sy})/sqrt(size(traces_vs{st,sy},1)),...
                nanmean(traces_vs{st,sy})-nanstd(traces_vs{st,sy})/sqrt(size(traces_vs{st,sy},1)),...
                col_resp(st), col_resp(st)); % Plot SEM
            hold on
            
            % Permutation test
            if st == numel(stnos)
                for cn = 1:size(compo,1)
                    c1 = compo(cn,1);
                    c2 = compo(cn,2);
                    sig_out = permTest_array(traces_vs{c1,sy}(~isnan(traces_vs{c1,sy}(:,1)),:),...
                        traces_vs{c2,sy}(~isnan(traces_vs{c2,sy}(:,1)),:),5000);
                    xw = xv{sy}; xw(sig_out>0.05)=nan;
                    yl = get(gca, 'ylim');
                    plot(xw, ones(numel(xw),1)*(yl(2)-0.03*cn), 'linewidth', 3, 'color', col_rep(cn))
                end
            end
            
            % Plot parameters
            if sy == 2 && st == numel(stnos)
                title([gnames{gr}])
                ylabel('dF/F')
                xlabel('Time')
                line([5 5], yl, 'color', 'k')
                line([10 10], yl, 'color', 'k')
                line([21.5 21.5], yl, 'color', 'k')
                if st == 11
                    line([12.5 12.5], yl, 'color','k', 'linestyle', lspec{2})
                    line([17.5 17.5], yl, 'color', 'k', 'linestyle',lspec{3})
                end
                set(gca, 'ylim', yl)
                legend(s, {'fixedITI' sname320{1}}, 'location', 'best', 'orientation', 'vertical')
            end
        end
    end
    supertitle([gnames{gr} newline])
end

% 3.22a) Compare signal parameters between vITI conditions within group
% 3.22b) Compare signal parameters between vSD conditions within group
comps = {'amplitude' 'dfiti' 'fstart' 'fcue' 'slope' 'slope1' 'slope2' ...
    'mean' 'mean1' 'mean2' 'halfweighttime' 'bias' 'preresp'}; % Calcium signal parameters
nf = nf+2;
traces_comps = cell(1,1);
oc = 1;
% ocspec = {1, 2};
T_OUT_3 = cell(1,1);
test_outputs = cell(1,1);

for gr = 1:2
    figure(nf+gr-1)
    for sy = 1:2
        for co = 1:numel(comps)
            pval = ones(1,2);
            % Compare early and late stage trial data for each parameter
            k=1;
            for ra = unique(meta_in(meta_in(:,1)==gr,2))'% Loop rats
                trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==9 & meta_in(:,5)==1); % Select trials
                trid2 = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==10 & meta_in(:,6)==1 & meta_in(:,5)==1); % Select trials
                
                try % By-pass empty comparisons (sometimes it depends on sync whether or not there is one)
                    traces_comps{sy,co}(k,1,gr)= ...
                        nanmean(calcium(sy).(comps{co})(trid,:)); % Store mean values for each stage
                    traces_comps{sy,co}(k,2,gr)= ...
                        nanmean(calcium(sy).(comps{co})(trid2,:)); % Store mean values for each stage
                catch
                    traces_comps{sy,co}(k,1,gr)=nan(1,1);
                    traces_comps{sy,co}(k,2,gr)=nan(1,1);
                end
                k=k+1;
            end
            
            traces_comps{sy,co}(traces_comps{sy,co}==0)=nan;
            
            subplot(ceil(sqrt(numel(comps))), round(sqrt(numel(comps))), co) % Subplot
            title(comps{co})
            plot(repmat((1:2)+3*(sy-1), size(traces_comps{sy,co},1), 1),traces_comps{sy,co}(:,:,gr),'.','color','k', 'markersize', 12) % Scatter plot of all values
            hold on
            plot((1:2)+3*(sy-1), nanmean(traces_comps{sy,co}(:,:,gr)), 'xk') % Plot mean
            plot(repmat((1:2)+3*(sy-1), size(traces_comps{sy,co},1),1)',traces_comps{sy,co}(:,:,gr)', 'color','k') % Plot individual rat traces
            set(gca, 'xlim', [0 6])
            title(comps{co})
  
            % Test
            [~,pval,~,s_out]=ttest(traces_comps{sy,co}(:,1,gr), traces_comps{sy,co}(:,2,gr));
            test_outputs{gr,sy}(:,co) = [round(s_out.tstat,4), s_out.df, round(pval,4)];
            yl = get(gca, 'ylim');
            if pval<0.05
                line([1 2]+3*(sy-1), [ yl(2)+(yl(2)-yl(1))*0.1  yl(2)+(yl(2)-yl(1))*0.1], 'color', 'k')
                text(1.5+3*(sy-1),  yl(2)+(yl(2)-yl(1))*0.1, '*','HorizontalAlignment', 'center')
                if pval<0.01
                    text(1.5+3*(sy-1),  yl(2)+(yl(2)-yl(1))*0.1, '**','HorizontalAlignment', 'center')
                    if pval<0.001
                        text(1.5+3*(sy-1),  yl(2)+(yl(2)-yl(1))*0.1, '***','HorizontalAlignment', 'center')
                    end
                end
            end
            if sy == 2
                set(gca, 'ylim', [yl(1), yl(2)+(yl(2)-yl(1))*0.2])
            end
            clear s_out
            supertitle([gnames{gr} newline])
            
        end
        T_OUT_3{gr,sy}=array2table(test_outputs{gr,sy}, 'VariableNames', comps);

    end
end


                
        
        
                % Compare signal with stage 9 (baselines and all conditions)
% Compare signal parameters between stage 9 and vSD/vITI
% Compare all with each other

% Compare fast and slow learners

% Compare post errors etc.

% Save figures
savefolder3 = ['D:\Data\FPCompiled_stages\Analysis\' exptname]; % folder where this analysis will be saved
figHandles = findobj('Type', 'figure'); % Gather all figure handles
set(figHandles, 'renderer', 'painters'); % Set renderer to 'painters' for proper further processing in Illustrator
set(figHandles, 'Position', get(0, 'Screensize')); % Set figures to full screen size (easier to process images in Illustrator)

for fign = 1:numel(figHandles)
    saveas(figHandles(fign), [savefolder3, '\Calcium-HighCog-',num2str(figHandles(fign).Number)], 'epsc'); % Save each figure as .EPSC file 
end


%% Script should work until here!
% Haven't rly tested stuff beyond. Work in progress

%% 3.x) Distribution of signal properties per rat, per stage
nf = nf+1;

figure(nf)
for gr = 1:2
    for st = 1:9
        for ra = unique(meta_in(meta_in(:,1)==gr,2))'
            for vn = 1:numel(vnames)
                trid = (meta_in(:,1)==gr & meta_in(:,2)==ra & meta_in(:,3)==st);
                subplot(numel(vnames),9,st+9*(vn-1))
                histogram(calcium(sy).(vnames{vn})(trid,:),35,'Normalization','pdf','DisplayStyle','stairs')
                hold on
            end
        end
    end
end



% Distribution of signal properties within stage, trained vs untrained

% Distribution of signal properties






%% Basic signal analysis
ncond = max(unique(meta_in(:,3)));
ratMeans = cell(2,ncond); ratErrorMeans = cell(2,ncond);
ratTrainedMeans = cell(2,ncond); ratTrainedErrorMeans = cell(2,ncond);
ratTrials=cell(2,11); ratErrorTrials=cell(2,11);
ratTrainedTrials=cell(2,11); ratTrainedErrorTrials=cell(2,11); 
dataset = data_in_z{1};
winsize = size(dataset,2);

% clf
for gr = 1:2
    nsize = unique(meta_in(meta_in(:,1)==gr,2)); % id of all rats in group
    for st = 1:ncond
        for ra = 1:numel(nsize)
            tridM = (meta_in(:,1)==gr & meta_in(:,2)==nsize(ra) & ...
                meta_in(:,3)==st & meta_in(:,5)==1 & meta_in(:,11)==0); % Correct untrained
            tridEM = (meta_in(:,1)==gr & meta_in(:,2)==nsize(ra) & ...
                meta_in(:,3)==st & meta_in(:,5)~=1 & meta_in(:,11)==0); % Error untrained
            tridTM = (meta_in(:,1)==gr & meta_in(:,2)==nsize(ra) & ...
                meta_in(:,3)==st & meta_in(:,5)==1 & meta_in(:,11)==1); % Correct trained
            tridTEM = (meta_in(:,1)==gr & meta_in(:,2)==nsize(ra) & ...
                meta_in(:,3)==st & meta_in(:,5)~=1 & meta_in(:,11)==1); % Error trained
            
            if sum(tridM)<5
                ratMeans{gr,st}(ra,:)=nan(1,winsize);
                data=nan(1,winsize);
            else
                data = dataset(tridM,:);
                ratMeans{gr,st}(ra,:)=nanmean(data);
            end
            
            if sum(tridEM)<5
                ratErrorMeans{gr,st}(ra,:)=nan(1,winsize);
                dataErrors=nan(1,winsize);
            else
                dataErrors = dataset(tridEM,:);
                ratErrorMeans{gr,st}(ra,:)=nanmean(dataErrors);
            end
            
            if sum(tridTM)<5
                ratTrainedMeans{gr,st}(ra,:)=nan(1,winsize);
                dataTrained=nan(1,winsize);
            else
                dataTrained = dataset(tridTM,:);
                ratTrainedMeans{gr,st}(ra,:)=nanmean(dataTrained);
            end
            
            if sum(tridTEM)<5
                ratTrainedErrorMeans{gr,st}(ra,:)=nan(1,winsize);
                dataTrainedErrors=nan(1,winsize);
            else
                dataTrainedErrors = dataset(tridTEM,:);
                ratTrainedErrorMeans{gr,st}(ra,:)=nanmean(dataTrainedErrors);
            end
            
            if isempty(ratTrials{gr,ra})
                ratTrials{gr,ra}=data;
                ratErrorTrials{gr,ra}=dataErrors;
                ratTrainedTrials{gr,ra}=dataTrained;
                ratTrainedErrorTrials{gr,ra}=dataTrainedErrors;
            else
                ratTrials{gr,ra}=[ratTrials{gr,ra}; nan(2,size(data,2)); data];
                ratErrorTrials{gr,ra}=[ratErrorTrials{gr,ra}; nan(2,size(dataErrors,2)); dataErrors];
                ratTrainedTrials{gr,ra}=[ratTrainedTrials{gr,ra}; nan(2,size(dataTrained,2)); dataTrained];
                ratTrainedErrorTrials{gr,ra}=[ratTrainedErrorTrials{gr,ra}; nan(2,size(dataTrainedErrors,2)); dataTrainedErrors];
            end
        end %rat
    end %stage
end %group

% Visualizations

% Does trained vs non trained signal vs errors differ across stages?
% Plot group means
%   group 1, stage 1, all compare --> 2 figures with 9 subplots?
for gr = 1:2
    figure(6+gr) % Figures 2 and 3
    for st = 1:9
        subplot(3,3,st)
        hold on
%         p1 = plot(movmean(nanmean(ratMeans{gr,st}),25), 'color', col_rep(1), 'linewidth',2);
% %         p2 = plot(movmean(nanmean(ratErrorMeans{gr,st}),25), 'color', col_rep(2), 'linewidth',2);
%         p3 = plot(movmean(nanmean(ratTrainedMeans{gr,st}),25), 'color', col_rep(3), 'linewidth',2);
%         p4 = plot(movmean(nanmean(ratTrainedErrorMeans{gr,st}),25), 'color', col_rep(4),'linewidth',1);
        p1 = plot(movmean(ratMeans{gr,st},2,25)', 'color', col_rep(1), 'linewidth',.5);
        p2 = plot(movmean(ratTrainedMeans{gr,st},2,25)', 'color', col_rep(3), 'linewidth',.5);
        xlim([1 15*15.89])
        yl = get(gca,'ylim');
        line([5*15.89 5*15.89], yl, 'color', 'k')
        line([10*15.89 10*15.89], yl, 'color', 'k')
    end %st
end %gr


% Plot stage means
figure(9) % Figure 4
for gr = 1:2
    for st = 1:9
        % Means of untrained trials
        subplot(3,2,gr)
        hold on
        plot(movmean(nanmean(ratMeans{gr,st}),25), 'color', [0.3 1-(st/9) st/9], 'linewidth',st/3);
        xlim([1 15*15.89])
        yl = get(gca,'ylim');
        line([5*15.89 5*15.89], yl, 'color', 'k')
        line([10*15.89 10*15.89], yl, 'color', 'k')

        % Means of trained trials
        subplot(3,2,gr+2)
        hold on
        plot(movmean(nanmean(ratTrainedMeans{gr,st}),25), 'color', [st/9 1-(st/9) .3], 'linewidth',st/3);
        xlim([1 15*15.89])
        yl = get(gca,'ylim');
        line([5*15.89 5*15.89], yl, 'color', 'k')
        line([10*15.89 10*15.89], yl, 'color', 'k')

        % Means of (trained-untrained) trials
        subplot(3,2,gr+4)
        hold on
        plot(movmean(nanmean(ratMeans{gr,st}-ratTrainedMeans{gr,st}),25), 'color', [0.3 1-(st/9) st/9], 'linewidth',st/3);
        xlim([1 15*15.89])
        yl = get(gca,'ylim');
        line([5*15.89 5*15.89], yl, 'color', 'k')
        line([10*15.89 10*15.89], yl, 'color', 'k')

    end %st
end %gr

%% Correlations between stages

wStart = ceil(5*15.89);
wEnd = ceil(12*15.89);
dataCorr = cell(1,1);

for gr=1:2
figure(4+gr)
    for ra=1:size(ratMeans{gr},1)
        for st=1:9
            dataCorr{gr,ra}(st,:) = ratMeans{gr,st}(ra,wStart:wEnd);
        end
        dcMean{gr}(:,ra) = nanmean(dataCorr{gr,ra},2);
        dcInt{gr}(:,ra) = trapz(dataCorr{gr,ra},2);
%         dataCorr{gr,ra}(isnan(dataCorr{gr,ra}))=0;
        dcCorr{gr}(:,:,ra)  = corrcoef(dataCorr{gr,ra}', 'rows','pairwise');
        [dcPeakA{gr}(:,ra),dcPeakT{gr}(:,ra)] =max(dataCorr{gr,ra},[],2);
    end
    subplot(2,2,gr)
    imagesc(nanmean(dcCorr{gr},3))
    caxis([0 0.75])
    
    subplot(2,2,3)
    hold on
    dcPeakA{gr}(dcPeakA{gr}==1)=nan;
    boxplot(dcPeakA{gr}','positions',(0.85:1:8.85)+(0.3*(gr-1)),'plotstyle','compact','colors',col_resp(gr));
    
    subplot(2,2,4)
    hold on
    dcPeakT{gr}(dcPeakT{gr}==1)=nan;
    boxplot(dcPeakT{gr}','positions',(0.85:1:8.85)+(0.3*(gr-1)),'plotstyle','compact','colors',col_resp(gr));
end

%% Variance during delay per stage
ncond = max(unique(meta_in(:,3)));
dataset = data_in{1};
winsize = size(dataset,2);
dataV = cell(1,1);
dataA = cell(1,1);
stageVar = cell(1,1);
stageVarD = cell(1,1);
sync_win = ceil(5*15.89):ceil(10*15.89);
base_win = ceil(1:4*15.89);
% clf

for gr=1:2 % Loop through groups
    nsize = unique(meta_in(meta_in(:,1)==gr,2)); % id of all rats in group
    for ra = 1:numel(nsize) % Loop through rats
        m=1;
        for st = 1:9 % Loop through stages
            tridV =  (meta_in(:,1)==gr & meta_in(:,2)==nsize(ra) & ...
                meta_in(:,3)==st);%  & meta_in(:,5)>=1 & meta_in(:,11)>=0);
            
            % Calculate variance during delay
            % Variance
            if sum(tridV)<5
                data = nan(1,winsize);
            else
                data = dataset(tridV,:);
            end
            
            % Individual trials
            dataV{gr}{ra,st} = var(data(:,sync_win),0,2); % Variance
            dataA{gr}{ra,st} = trapz(data(:,sync_win),2); % AUC
            
            % Information distance
            binsize = ceil(numel(sync_win)/100);
            binnum = ceil(1:numel(sync_win)/max(binsize));
            % fraction of bins
            binfrac = 1/max(binnum);
            H = []; binmean =[];
            for trial = 1:size(data,1)
                if sync == 1
                    binmean(trial,:) = accumarray(binnum(:), ...
                        data(trial,sync_win)',...
                        [],@nanmean);
                    % overall mean signal
                    %                         trMed = nanmean(plotDat(trial,:));
                    %                 trMed = stageData.(allstages{stage}).baseWin(trial,3);
                else
                    binmean(trial,:) = accumarray(binnum(:), ...
                        data(trial,sync_win)',...
                        [],@nanmean);
                    %                 trMed = nanmean(plotDat(trial,:));
                    %                         normalized with baseline w
                    %                 trMed = stageData.(allstages{stage}).baseWin(trial,3);
                end
                
                for bin = 1:max(binnum)
                    %             H(trial, bin)= binfrac*binmean(trial,bin)*(log2(binmean(trial,bin)/tmp1(trial)));
                    H(trial, bin)= binfrac*binmean(trial,bin)*(log2(abs(binmean(trial,bin)-nanmedian(data(trial,base_win)))));
                end
                
            end
            
            
        end
    end
end

%             
%             
%             figure(6+gr)
%             subplot(floor(sqrt(numel(nsize))),ceil(sqrt(numel(nsize))),ra)
%             hold on;
%             plotdat = dataA{gr}{ra,st};
%             xF = m:m-1+size(plotdat,1);
%             plot(xF, plotdat,'.','color',col_rep(st))
%             coeffs = polyfit(1:size(plotdat,1),plotdat,1);
%             yF = polyval(coeffs,1:size(plotdat,1));
%             stageVar{gr}(ra,st) = mean(plotdat);
%             stageVarD{gr}(ra,st) = coeffs(1);
%             plot(xF,yF,'color','k')
%             m=m+size(plotdat,1);
%         end
%     end
%     figure(9)
%     subplot(2,2,gr)
%     plot(repmat(1:9,numel(nsize),1),stageVar{gr},'.', 'markersize',12)
%     xlim([0 10])
%     
%     subplot(2,2,gr+2)
%     plot(repmat(1:9,numel(nsize),1),stageVarD{gr},'.', 'markersize',12)
%     xlim([0 10])
% end

