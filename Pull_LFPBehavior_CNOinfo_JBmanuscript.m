% Function to pull and save LFP and behavior (statetimes) info from
% requested animals for specified time periods after CNO injections and
% during baseline periods (for further analysis later)
% Subfunction of Analyze_LFPBehavior_CNO_JBmanuscript.m

% Pulls data from 12 hour periods after specified CNO times and baseline
% times, so that can analyze variable smaller length periods later

% INPUT: 
%       - AnimNames is a cell array of animal names (as character strings) you want to pull

% OUTPUT: 
%       - LFPBehavior_Animals is a structure that contains each input animal's
%       saved structure within it
%       - Saves each animal's individual info (so don't have to load huge
%       LFP file again) as AnimName_LFPBehavior_CNOinfo.mat 

%State codes:
% 1 = REM
% 2 = NREM
% 3 = Gen Wake ** NOT USED **
% 4 = Active Wake
% 5 = Quiet Wake

function LFPBehavior_Animals = Pull_LFPBehavior_CNOinfo_JBmanuscript(AnimNames)


LFPBehavior_Animals = struct;

for i = 1:length(AnimNames)

    % load LFPinfo and statetimes for each animal
    thisAnim = char(AnimNames(i));
    LFPfile = dir([thisAnim, '_LFPinfo*']);
    LFPinfo = load(LFPfile.name);    
    LFPinfo = LFPinfo.LFPinfo;
    statetimes = load([thisAnim, '_STATETIMES.mat']);

    if isstruct(statetimes)
        if isfield(statetimes, 'statetimes')
            statetimes = statetimes.statetimes;
        else
            statetimes = statetimes.STATETIMES;
        end
    end

    
    fprintf('Starting LFP prep for %s... \n', thisAnim)
    % PREP LFPinfo, pull relevant data, frequency bands, align times to
    % behavior statetimes
    
    DeltaBand = [.5 4];
    ThetaBand = [6 8];
    AlphaBand = [9 15];
    
    num_blocks = size (LFPinfo);
    num_blocks = num_blocks(2);
    % skip any short blocks at the beginning
    for jj=1:num_blocks
        if LFPinfo(jj).duration>3599
            disp ('start block:')
            disp (jj)
            break
        else
            continue
        end
    end
    
    % find start time of LFP data to align with behavior
    start_time_string = datestr(unixtime(LFPinfo(jj).startTime));
    LFPhour=str2num(start_time_string(end-7:end-6));
    LFPminute=str2num(start_time_string(end-4:end-3));
    LFPseconds=str2num(start_time_string(end-1:end));

    time_from_midnight=(LFPhour*3600)+(LFPminute*60)+LFPseconds; % in seconds: time from beginning of BL1: start of first block
    LFP_start = time_from_midnight - (7.5*3600); % time expt started in seconds from 7.30 am that first day (BL1)

    % pull out delta, theta, and alpha power of full experiment
    total_times = [];
    total_delta_powers = []; % SAMPLING RATE: TWO POINTS EVERY SECOND!!
    total_theta_powers = [];
    total_alpha_powers = [];
    TOTAL_powers = [];

    myStart = jj;
    if thisAnim == 'JB13' | thisAnim == 'JB14'
        myEnd = 133;
    else
        myEnd = 200;
    end
    % loop through each hour block at a time
    for ii=myStart:myEnd

    total_powers = sum(LFPinfo(ii).spectrogram.P);
    TOTAL_powers = [TOTAL_powers, total_powers];
    lower_delta = find (LFPinfo(ii).spectrogram.F>=DeltaBand(1), 1); 
    upper_delta = find (LFPinfo(ii).spectrogram.F>=DeltaBand(2), 1);

    lower_theta = find(LFPinfo(ii).spectrogram.F>=ThetaBand(1), 1);
    upper_theta = find(LFPinfo(ii).spectrogram.F>=ThetaBand(2), 1);

    delta_power = LFPinfo(ii).spectrogram.P(lower_delta:upper_delta,:);
    delta_power = sum(delta_power);
    delta_power = delta_power./total_powers; % pulling out relative bands here, normalized to total power at that time 

    theta_power = LFPinfo(ii).spectrogram.P(lower_theta:upper_theta,:);
    theta_power = sum(theta_power);
    theta_power = theta_power./total_powers;

    lower_alpha = find(LFPinfo(ii).spectrogram.F>=AlphaBand(1), 1);
    upper_alpha = find(LFPinfo(ii).spectrogram.F>=AlphaBand(2), 1);

    alpha_power = LFPinfo(ii).spectrogram.P(lower_alpha:upper_alpha,:);
    alpha_power = sum(alpha_power);
    alpha_power = alpha_power./total_powers;

    % pull out timestamps of each LFP point
    this_start = ((ii-myStart)*3600)+ LFP_start; 
    theseTimes = [(LFPinfo(ii).spectrogram.T) + this_start];
    total_times = [total_times, theseTimes];

    total_delta_powers = [total_delta_powers, delta_power];
    total_theta_powers = [total_theta_powers, theta_power];
    total_alpha_powers = [total_alpha_powers, alpha_power];
    end

    %%% Convert statetimes from unix time to seconds from 7:30 am on day
    %%% 1 of recording

    unixtimes = statetimes(:,2);
    mytimes = unixtime(unixtimes);
    day_start = mytimes(1, 3);
    days = mytimes(:, 3);
    days = days - day_start;
    hours = mytimes(:,4);
    minutes = mytimes(:,5);
    seconds = mytimes(:,6);
    totalSeconds = (days*24*3600) + (hours*3600) + (minutes*60) + seconds;
    MYseconds = totalSeconds - (7.5*3600); % times in seconds since 7:30 am on first day of recording
    % if experiment spans 2 months
    if length(unique(mytimes(:,2)))>1
        months = mytimes(:,2);
        days = mytimes(:,3);

        if ismember(31, days)
            x = 31;
        elseif ismember(30, days)
            x = 30;
        elseif ismember(29, days)
            x = 29;
        elseif ismember(28, days)
            x = 28;
        end

        for d = 1:days(end)
            days(find(days==d)) = x+d;
        end
        days = days - day_start;
        totalSeconds = (days*24*3600) + (hours*3600) + (minutes*60) + seconds;
        MYseconds = totalSeconds - (7.5*3600); % times in seconds since 7:30 am on first day of recording
    end

    new_statetimes = statetimes(:,1);
    new_statetimes(:,2) = MYseconds;

    % get rid of duplicate states
     testD   = diff(new_statetimes(:,1));
     killem  = find(testD == 0);
     killem  = killem+1;
     new_statetimes(killem,:) = [];

    % find and eliminate transient wake episodes (e.g. the animal wakes up
    % simply to shift position)
    nsec=60; %set limit at 1 min-eliminate all states that only last 1 minute or less
    sdur = diff(new_statetimes(:,2));
    fcks = find(new_statetimes(1:end-1,1)>2 & sdur<=nsec ); % find instances of waking that are less than nsec seconds long

    tempf = new_statetimes(fcks+1,1)<3; % find one fcks for which the following state is still sleeping
    fcks2 = fcks(tempf);
    new_statetimes(fcks2,:) = [];
    % second pass
    testD   = diff(new_statetimes(:,1));
    killem  = find(testD == 0);
    killem  = killem+1;
    new_statetimes(killem,:) = [];
    
    %%%%%%%%%%%%%% Now start Pulling relevant data to save! %%%%%%%%%%%%%%%%%%%

    %%%% DEFINE LENGTH of time to pull - will refine later,
    %%%% pulling max period now
    CNO_period = 12; % in hours
    CNO_period = CNO_period*3600; % in seconds

    %%%%% ENTER CNO TIMES FOR EACH ANIMAL %%%%%%%%

    if thisAnim == 'JB13' | thisAnim == 'JB14'
        CNO_start = ((24*3) + 1.75) * 3600;
        CNO = [CNO_start, CNO_start + 10.25*3600, CNO_start + 23*3600, CNO_start + 34.25*3600];
        CNO_light =[CNO(1), CNO(3)];
        CNO_dark = [CNO(2), CNO(4)];
    elseif thisAnim == 'JB22' | thisAnim == 'JB21'
        CNO_start = 24*5*3600;
        CNO = [CNO_start, CNO_start + (12*3600), CNO_start + (24*3600), CNO_start + (36*3600)];
        CNO_light = [CNO(1)]; % MD4 light period neural signal cut out, so don't trust LFP
        CNO_dark = [CNO(2)];
    elseif thisAnim == 'JB26' | thisAnim == 'JB27' | thisAnim == 'JB32' | thisAnim == 'JB33'
        CNO_start = 24*5*3600;
        CNO = [CNO_start, CNO_start + (12*3600), CNO_start + (24*3600), CNO_start + (36*3600)];
        CNO_light = [CNO(1), CNO(3)]; 
        CNO_dark = [CNO(2), CNO(4)];
    else
        disp ('ERROR! ENTER CNO INJECTION TIMES FOR THIS ANIMAL AND TRY AGAIN!')
        disp (thisAnim)
        keyboard
    end
    
    %%%%%%%%% ENTER BASELINE TIMES  %%%%%%%%%%%%

    if thisAnim == 'JB13' 
        baseline_lightperiods_start = [ 24*2*3600, 24*5*3600];
        baseline_darkperiods_start = [2.5*24*3600, 5.5*24*3600];      
    elseif thisAnim == 'JB14'
        % NOTE: this animal had messed up signal @ day 1-1.2ish and 5-5.2ish. Be careful with Light Period baseline data
        baseline_lightperiods_start = [24*2*3600, 24*5.2*3600];
        baseline_darkperiods_start = [2.5*24*3600, 5.5*24*3600];  
    elseif thisAnim == 'JB21' | thisAnim == 'JB22'
        baseline_lightperiods_start = [24*3*3600, 24*4*3600]; % only one period to match only one CNO period being analyzed
        baseline_darkperiods_start = [3.5*24*3600, 4.5*24*3600]; % only one period to match only one CNO period being analyzed
    elseif thisAnim == 'JB26' | thisAnim == 'JB27' | thisAnim == 'JB32' | thisAnim == 'JB33'
        baseline_lightperiods_start = [24*4*3600, 24*3*3600];
        baseline_darkperiods_start = [4.5*24*3600, 3.5*24*3600];
    else
        disp ('ERROR! ENTER CNO INJECTION TIMES FOR THIS ANIMAL AND TRY AGAIN!')
        disp (thisAnim)
        keyboard
    end

    % Go Pull out relevant LFP and statetimes for each of these periods,
    % save in this animal's structure
    num_periods = length(CNO_light) + length(CNO_dark) + length(baseline_lightperiods_start) + length(baseline_darkperiods_start);
    clear LFPBehavior
    LFPBehavior = struct;
    
    for uu = 1:num_periods

        if uu <= length(CNO_light)
            period = [CNO_light(uu), CNO_light(uu)+CNO_period];
            these_statetimes = new_statetimes(find(new_statetimes(:,2)>=period(1)),:);
            these_statetimes = these_statetimes(find(these_statetimes(:,2)<=period(2)),:);
            first_statetime = new_statetimes(find(new_statetimes(:,2)<=period(1), 1, 'last'), 1);
            these_statetimes = [first_statetime, period(1); these_statetimes];
            LFPBehavior(uu).PeriodType = {'CNO Light'};

        elseif uu <= (length(CNO_dark) + length(CNO_light))
            period = [CNO_dark(uu-length(CNO_light)), CNO_dark(uu-length(CNO_light))+CNO_period];
            these_statetimes = new_statetimes(find(new_statetimes(:,2)>=period(1)),:);
            these_statetimes = these_statetimes(find(these_statetimes(:,2)<=period(2)),:);
            first_statetime = new_statetimes(find(new_statetimes(:,2)<=period(1), 1, 'last'), 1);
            these_statetimes = [first_statetime, period(1); these_statetimes];
            LFPBehavior(uu).PeriodType = {'CNO Dark'};

        elseif uu <= (length(baseline_lightperiods_start) + length(CNO_dark) + length(CNO_light))
            period = [baseline_lightperiods_start(uu-length(CNO_light) - length(CNO_dark)), baseline_lightperiods_start(uu-length(CNO_light) - length(CNO_dark))+CNO_period];
            these_statetimes = new_statetimes(find(new_statetimes(:,2)>=period(1)),:);
            these_statetimes = these_statetimes(find(these_statetimes(:,2)<=period(2)),:);       
            first_statetime = new_statetimes(find(new_statetimes(:,2)<=period(1), 1, 'last'), 1);
            these_statetimes = [first_statetime, period(1); these_statetimes];
            LFPBehavior(uu).PeriodType = {'Baseline Light'};

        elseif uu <= (length(baseline_darkperiods_start) + length(baseline_lightperiods_start) + length(CNO_dark) + length(CNO_light))
            period = [baseline_darkperiods_start(uu-length(CNO_light)-length(CNO_dark)-length(baseline_lightperiods_start)), baseline_darkperiods_start(uu-length(CNO_light)-length(CNO_dark)-length(baseline_lightperiods_start))+CNO_period];
            these_statetimes = new_statetimes(find(new_statetimes(:,2)>=period(1)),:);
            these_statetimes = these_statetimes(find(these_statetimes(:,2)<=period(2)),:);
            first_statetime = new_statetimes(find(new_statetimes(:,2)<=period(1), 1, 'last'), 1);
            these_statetimes = [first_statetime, period(1); these_statetimes];
            LFPBehavior(uu).PeriodType = {'Baseline Dark'};
        end

        % starting the period at the first behavioral transition after the start
        % of each period - the first behavioral state at the time of the
        % injections is probably artifact anyway (eg. stressful injection)
        total_time = period(2)-these_statetimes(1,2);

        freqstart = find (total_times>=these_statetimes(1,2), 1);
        freqstop = find(total_times<=period(2),1, 'last');
        total_delta = total_delta_powers(freqstart:freqstop);
        total_theta = total_theta_powers(freqstart:freqstop);
        total_alpha = total_alpha_powers(freqstart:freqstop);
        
        LFPBehavior(uu).total_times = total_times(freqstart:freqstop);
        LFPBehavior(uu).total_delta = total_delta;
        LFPBehavior(uu).total_theta = total_theta;
        LFPBehavior(uu).total_alpha = total_alpha;
        LFPBehavior(uu).statetimes = these_statetimes;
        LFPBehavior(uu).Period = period;

    end
    
    LFPBehavior_Animals.(thisAnim) = LFPBehavior;


end

end
