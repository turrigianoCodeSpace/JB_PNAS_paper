% Function to analyze LFPBehavior data for requested animals over specified
% amount of time after CNO injections compared to baseline periods
% For ACh manuscript Figure 1F,G; Table S1

% Analyzes specified length periods after CNO injections and baseline (time
% matched) periods from previously pulled 12 hour periods (if animal hasn't
% been analyzed yet, it calls the function to pull the data for that animal and saves it)
% Also plots combined data by animal and by period (calls separate function
% for the plotting, using only info in the output LFPBehavior_CombinedData
% structure)

% INPUT: 
%       - Animals is a cell array of animal names (as character strings)
%       you want to anlayze
%       - Period is the length of time to analyze after each CNO injection
%       and each baseline period (in hours)
%       - Optional third variable is START of the period to analyze (in hours, after start of each period, defined as time of CNO injection and time matched baseline period) - this
%       is for making time courses, so can bin in shorter chunks and
%       analyze different chunks at a time. If no input, default is 0
%       (analysis starts at the beginning of each period
%       - Optional 4th variable is whether or not to feed this time period
%       analysis into Plot_LFPBehavior_CNO_LightDarkCombine_JBmanuscript.m, which plots CNO data from
%       this one period normalized to time matched baseline period.

% OUTPUT: 
%       - LFPBehavior_CombinedData is a structure used to plot analyzed
%       data from all requested animals, with data organized by animal and
%       by period. This structure can also be used to generate stats
%       - Can save this structure if you want, but doesn't do it
%       automatically - has field for length of periods being analyzed

% PLOTS: 
%       - Behavioral data (percent time, epoch length etc for REM, NREM, AW, QW) after CNO,
%       normalized to baseline
%       - LFP data (delta, theta, alpha power in each behavioral state)
%       after CNO, normalized to baseline period
%       - Plots these things by animal and by CNO period

% input for ACh manuscript: ({'JB13', 'JB14', 'JB22', 'JB26', 'JB32', 'JB33', 8})
function [LFPBehavior_CombinedData, stats] = Analyze_LFPBehavior_CNO_LightDarkCombine_JBmanuscript(Animals, Period, varargin)

if nargin == 2
    StartPeriod = 0;
    pplot = 1;
elseif nargin == 3
    StartPeriod = varargin{1}*3600;
    pplot = 1;
else
    StartPeriod = varargin{1}*3600;
    pplot = varargin{2};
end

LFPBehavior_CombinedData = struct;
LFPBehavior_CombinedData.PeriodLength = Period;
LFPBehavior_CombinedData.StartPeriod = StartPeriod;
all_Animals = struct;

% pull each individual animal's LFP and behavior data

    [more_Animals] = Pull_LFPBehavior_CNOinfo_JBmanuscript(Animals);
    for n = 1:length(Animals)
        all_Animals.(char(Animals(n))).LFPBehavior = more_Animals.(char(Animals(n)));
    end
    
fprintf('Pulled all animal LFPBehavior data, Starting Analysis \n')

% Loop through all animals 
for AA = 1:length(Animals)
    thisAnimal = char(Animals(AA));
    fprintf('Starting %s... \n', thisAnimal)
    
    % Loop through each period (CNO and baseline periods)
    % keep track of how many of each type of period there are for the 'by
    % animal' analysis 
    CNO_light_num = 0;
    CNO_dark_num = 0;
    BL_light_num = 0;
    BL_dark_num = 0;
    for PP = 1:length(all_Animals.(thisAnimal).LFPBehavior)
        period = all_Animals.(thisAnimal).LFPBehavior(PP).Period;
        period(1) = period(1) + StartPeriod;
        period(2) = period(1)+(Period*3600); % new period length to actually analyze
        these_statetimes = all_Animals.(thisAnimal).LFPBehavior(PP).statetimes;
        these_statetimes = these_statetimes(find(these_statetimes(:,2)>=period(1)),:);
        these_statetimes = these_statetimes(find(these_statetimes(:,2)<=period(2)),:);
        
        total_times = all_Animals.(thisAnimal).LFPBehavior(PP).total_times;
        total_delta_powers = all_Animals.(thisAnimal).LFPBehavior(PP).total_delta;
        total_theta_powers = all_Animals.(thisAnimal).LFPBehavior(PP).total_theta;
        total_alpha_powers = all_Animals.(thisAnimal).LFPBehavior(PP).total_alpha;
        Period_type = all_Animals.(thisAnimal).LFPBehavior(PP).PeriodType;
        if strcmp(char(Period_type), 'CNO Light')
            CNO_light_num = CNO_light_num + 1;
        elseif strcmp(char(Period_type), 'CNO Dark')
            CNO_dark_num = CNO_dark_num + 1;
        elseif strcmp(char(Period_type), 'Baseline Light')
            BL_light_num = BL_light_num + 1;
        elseif strcmp(char(Period_type), 'Baseline Dark')
            BL_dark_num = BL_dark_num  + 1;
        end
        
        total_time = period(2)-period(1);
        
        freqstart = find(total_times>=period(1),1);
        freqstop = find(total_times<=period(2),1, 'last');
        total_delta = mean(total_delta_powers(freqstart:freqstop));
        total_theta = mean(total_theta_powers(freqstart:freqstop));
        total_alpha = mean(total_alpha_powers(freqstart:freqstop));

            REM = find (these_statetimes(:,1)==1);
            NREM = find (these_statetimes(:,1)==2);
            AW = find (these_statetimes(:,1)==4);
            QW = find (these_statetimes(:,1)==5);
            
        % length of each epoch
        allREMtimes = [];
        allNREMtimes = [];
        allAWtimes = [];
        allQWtimes = [];

        % these are average delta/theta/alpha power in each epoch
        delta_REM = [];
        theta_REM = [];
        alpha_REM = [];

        delta_NREM = [];
        theta_NREM = [];
        alpha_NREM = [];

        delta_AW = [];
        theta_AW = [];
        alpha_AW = [];

        delta_QW = [];
        theta_QW = [];
        alpha_QW = [];
    
        for r = 1:length(REM)

            if REM(r) == length(these_statetimes(:,1))
                thisEpoch = [these_statetimes(REM(r),2), period(2)]; 
            else
                thisEpoch = [these_statetimes(REM(r),2), these_statetimes(REM(r)+1,2)];                
            end
                start_freqanalysis = find(total_times>=thisEpoch(1), 1);
                stop_freqanalysis = find(total_times<=thisEpoch(2),1, 'last');

                delta_REM = [delta_REM, mean(total_delta_powers(start_freqanalysis:stop_freqanalysis))];
                theta_REM = [theta_REM, mean(total_theta_powers(start_freqanalysis:stop_freqanalysis))];
                alpha_REM = [alpha_REM, mean(total_alpha_powers(start_freqanalysis:stop_freqanalysis))];

            thisLength = thisEpoch(2)-thisEpoch(1);
            allREMtimes = [allREMtimes, thisLength];
        end

        for n = 1:length(NREM)

            if NREM(n) == length(these_statetimes(:,1))
                thisEpoch = [these_statetimes(NREM(n),2), period(2)];
            else 
                thisEpoch = [these_statetimes(NREM(n),2), these_statetimes(NREM(n)+1,2)];
            end
                start_freqanalysis = find (total_times>=thisEpoch(1), 1);
                stop_freqanalysis = find(total_times<=thisEpoch(2),1, 'last');

                delta_NREM = [delta_NREM, mean(total_delta_powers(start_freqanalysis:stop_freqanalysis))];
                theta_NREM = [theta_NREM, mean(total_theta_powers(start_freqanalysis:stop_freqanalysis))];
                alpha_NREM = [alpha_NREM, mean(total_alpha_powers(start_freqanalysis:stop_freqanalysis))];

            thisLength = thisEpoch(2)-thisEpoch(1);
            allNREMtimes = [allNREMtimes, thisLength];
        end

        for a = 1:length(AW)

            if AW(a) == length(these_statetimes(:,1))
                thisEpoch = [these_statetimes(AW(a),2), period(2)];
            else 
                thisEpoch = [these_statetimes(AW(a),2), these_statetimes(AW(a)+1,2)];
            end
                start_freqanalysis = find (total_times>=thisEpoch(1), 1);
                stop_freqanalysis = find(total_times<=thisEpoch(2),1, 'last');

                delta_AW = [delta_AW, mean(total_delta_powers(start_freqanalysis:stop_freqanalysis))];
                theta_AW = [theta_AW, mean(total_theta_powers(start_freqanalysis:stop_freqanalysis))];
                alpha_AW = [alpha_AW, mean(total_alpha_powers(start_freqanalysis:stop_freqanalysis))];

            thisLength = thisEpoch(2)-thisEpoch(1);
            allAWtimes = [allAWtimes, thisLength];
        end

        for q = 1:length(QW)

            if QW(q) == length(these_statetimes(:,1))
                thisEpoch = [these_statetimes(QW(q),2), period(2)];
            else 
                thisEpoch = [these_statetimes(QW(q),2), these_statetimes(QW(q)+1,2)];
            end
                start_freqanalysis = find (total_times>=thisEpoch(1), 1);
                stop_freqanalysis = find(total_times<=thisEpoch(2),1, 'last');

                delta_QW = [delta_QW, mean(total_delta_powers(start_freqanalysis:stop_freqanalysis))];
                theta_QW = [theta_QW, mean(total_theta_powers(start_freqanalysis:stop_freqanalysis))];
                alpha_QW = [alpha_QW, mean(total_alpha_powers(start_freqanalysis:stop_freqanalysis))];

            thisLength = thisEpoch(2)-thisEpoch(1);
            allQWtimes = [allQWtimes, thisLength];
        end

        REMtime = sum(allREMtimes);
        NREMtime = sum(allNREMtimes);
        AWtime = sum(allAWtimes);
        QWtime = sum(allQWtimes);
    
    % store behavior and LFP info for each period 
    if strcmp(char(Period_type), 'CNO Light')
        if ~isfield(LFPBehavior_CombinedData, 'DeltaPeriod_light_CNO')
            LFPBehavior_CombinedData.DeltaPeriod_light_CNO = [nanmean(delta_REM), nanmean(delta_NREM), nanmean(delta_AW), nanmean(delta_QW)];
            LFPBehavior_CombinedData.ThetaPeriod_light_CNO = [nanmean(theta_REM), nanmean(theta_NREM), nanmean(theta_AW), nanmean(theta_QW)];
            LFPBehavior_CombinedData.AlphaPeriod_light_CNO = [nanmean(alpha_REM), nanmean(alpha_NREM), nanmean(alpha_AW), nanmean(alpha_QW)];
            LFPBehavior_CombinedData.totalsPeriod_light_CNO = [total_delta, total_theta, total_alpha];
            LFPBehavior_CombinedData.PercentsPeriod_light_CNO = [REMtime/total_time, NREMtime/total_time, AWtime/total_time, QWtime/total_time];
            LFPBehavior_CombinedData.DurationsPeriod_light_CNO = [nanmean(allREMtimes), nanmean(allNREMtimes), nanmean(allAWtimes), nanmean(allQWtimes)];
            
            LFPBehavior_CombinedData.DeltaAnimal_light_CNO = [];
            LFPBehavior_CombinedData.ThetaAnimal_light_CNO = [];
            LFPBehavior_CombinedData.AlphaAnimal_light_CNO = [];
            LFPBehavior_CombinedData.totalsAnimal_light_CNO = [];
            LFPBehavior_CombinedData.PercentsAnimal_light_CNO = [];
            LFPBehavior_CombinedData.DurationsAnimal_light_CNO = [];
            LFPBehavior_CombinedData.ANIMALS_CNOlight = {};
        else
            LFPBehavior_CombinedData.DeltaPeriod_light_CNO(end+1, :) = [nanmean(delta_REM), nanmean(delta_NREM), nanmean(delta_AW), nanmean(delta_QW)];
            LFPBehavior_CombinedData.ThetaPeriod_light_CNO(end+1, :) = [nanmean(theta_REM), nanmean(theta_NREM), nanmean(theta_AW), nanmean(theta_QW)];
            LFPBehavior_CombinedData.AlphaPeriod_light_CNO(end+1, :) = [nanmean(alpha_REM), nanmean(alpha_NREM), nanmean(alpha_AW), nanmean(alpha_QW)];
            LFPBehavior_CombinedData.totalsPeriod_light_CNO(end+1, :) = [total_delta, total_theta, total_alpha];
            LFPBehavior_CombinedData.PercentsPeriod_light_CNO(end+1, :) = [REMtime/total_time, NREMtime/total_time, AWtime/total_time, QWtime/total_time];
            LFPBehavior_CombinedData.DurationsPeriod_light_CNO(end+1,:) = [nanmean(allREMtimes), nanmean(allNREMtimes), nanmean(allAWtimes), nanmean(allQWtimes)];
        end
  
    elseif strcmp(char(Period_type), 'CNO Dark')
        if ~isfield(LFPBehavior_CombinedData, 'DeltaPeriod_dark_CNO')
            LFPBehavior_CombinedData.DeltaPeriod_dark_CNO = [nanmean(delta_REM), nanmean(delta_NREM), nanmean(delta_AW), nanmean(delta_QW)];
            LFPBehavior_CombinedData.ThetaPeriod_dark_CNO = [nanmean(theta_REM), nanmean(theta_NREM), nanmean(theta_AW), nanmean(theta_QW)];
            LFPBehavior_CombinedData.AlphaPeriod_dark_CNO = [nanmean(alpha_REM), nanmean(alpha_NREM), nanmean(alpha_AW), nanmean(alpha_QW)];
            LFPBehavior_CombinedData.totalsPeriod_dark_CNO = [total_delta, total_theta, total_alpha];
            LFPBehavior_CombinedData.PercentsPeriod_dark_CNO = [REMtime/total_time, NREMtime/total_time, AWtime/total_time, QWtime/total_time];
            LFPBehavior_CombinedData.DurationsPeriod_dark_CNO = [nanmean(allREMtimes), nanmean(allNREMtimes), nanmean(allAWtimes), nanmean(allQWtimes)];
                    
            LFPBehavior_CombinedData.DeltaAnimal_dark_CNO = [];
            LFPBehavior_CombinedData.ThetaAnimal_dark_CNO = [];
            LFPBehavior_CombinedData.AlphaAnimal_dark_CNO = [];
            LFPBehavior_CombinedData.totalsAnimal_dark_CNO = [];
            LFPBehavior_CombinedData.PercentsAnimal_dark_CNO = [];
            LFPBehavior_CombinedData.DurationsAnimal_dark_CNO = [];
            LFPBehavior_CombinedData.ANIMALS_CNOdark = {};
        else
            LFPBehavior_CombinedData.DeltaPeriod_dark_CNO(end+1, :) = [nanmean(delta_REM), nanmean(delta_NREM), nanmean(delta_AW), nanmean(delta_QW)];
            LFPBehavior_CombinedData.ThetaPeriod_dark_CNO(end+1, :) = [nanmean(theta_REM), nanmean(theta_NREM), nanmean(theta_AW), nanmean(theta_QW)];
            LFPBehavior_CombinedData.AlphaPeriod_dark_CNO(end+1, :) = [nanmean(alpha_REM), nanmean(alpha_NREM), nanmean(alpha_AW), nanmean(alpha_QW)];
            LFPBehavior_CombinedData.totalsPeriod_dark_CNO(end+1, :) = [total_delta, total_theta, total_alpha];
            LFPBehavior_CombinedData.PercentsPeriod_dark_CNO(end+1, :) = [REMtime/total_time, NREMtime/total_time, AWtime/total_time, QWtime/total_time];
            LFPBehavior_CombinedData.DurationsPeriod_dark_CNO(end+1,:) = [nanmean(allREMtimes), nanmean(allNREMtimes), nanmean(allAWtimes), nanmean(allQWtimes)];
        end
      
    elseif strcmp(char(Period_type), 'Baseline Light')
        if ~isfield(LFPBehavior_CombinedData, 'DeltaPeriod_light_baseline')
            LFPBehavior_CombinedData.DeltaPeriod_light_baseline = [nanmean(delta_REM), nanmean(delta_NREM), nanmean(delta_AW), nanmean(delta_QW)];
            LFPBehavior_CombinedData.ThetaPeriod_light_baseline = [nanmean(theta_REM), nanmean(theta_NREM), nanmean(theta_AW), nanmean(theta_QW)];
            LFPBehavior_CombinedData.AlphaPeriod_light_baseline = [nanmean(alpha_REM), nanmean(alpha_NREM), nanmean(alpha_AW), nanmean(alpha_QW)];
            LFPBehavior_CombinedData.totalsPeriod_light_baseline = [total_delta, total_theta, total_alpha];
            LFPBehavior_CombinedData.PercentsPeriod_light_baseline = [REMtime/total_time, NREMtime/total_time, AWtime/total_time, QWtime/total_time];
            LFPBehavior_CombinedData.DurationsPeriod_light_baseline = [nanmean(allREMtimes), nanmean(allNREMtimes), nanmean(allAWtimes), nanmean(allQWtimes)];
            
            LFPBehavior_CombinedData.DeltaAnimal_light_baseline = [];
            LFPBehavior_CombinedData.ThetaAnimal_light_baseline = [];
            LFPBehavior_CombinedData.AlphaAnimal_light_baseline = [];
            LFPBehavior_CombinedData.totalsAnimal_light_baseline = [];
            LFPBehavior_CombinedData.PercentsAnimal_light_baseline = [];
            LFPBehavior_CombinedData.DurationsAnimal_light_baseline = [];
            LFPBehavior_CombinedData.ANIMALS_BLlight = {};
        else
            LFPBehavior_CombinedData.DeltaPeriod_light_baseline(end+1, :) = [nanmean(delta_REM), nanmean(delta_NREM), nanmean(delta_AW), nanmean(delta_QW)];
            LFPBehavior_CombinedData.ThetaPeriod_light_baseline(end+1, :) = [nanmean(theta_REM), nanmean(theta_NREM), nanmean(theta_AW), nanmean(theta_QW)];
            LFPBehavior_CombinedData.AlphaPeriod_light_baseline(end+1, :) = [nanmean(alpha_REM), nanmean(alpha_NREM), nanmean(alpha_AW), nanmean(alpha_QW)];
            LFPBehavior_CombinedData.totalsPeriod_light_baseline(end+1, :) = [total_delta, total_theta, total_alpha];
            LFPBehavior_CombinedData.PercentsPeriod_light_baseline(end+1, :) = [REMtime/total_time, NREMtime/total_time, AWtime/total_time, QWtime/total_time];
            LFPBehavior_CombinedData.DurationsPeriod_light_baseline(end+1,:) = [nanmean(allREMtimes), nanmean(allNREMtimes), nanmean(allAWtimes), nanmean(allQWtimes)];
        end
     
    elseif strcmp(char(Period_type), 'Baseline Dark')
        if ~isfield(LFPBehavior_CombinedData, 'DeltaPeriod_dark_baseline')
            LFPBehavior_CombinedData.DeltaPeriod_dark_baseline = [nanmean(delta_REM), nanmean(delta_NREM), nanmean(delta_AW), nanmean(delta_QW)];
            LFPBehavior_CombinedData.ThetaPeriod_dark_baseline = [nanmean(theta_REM), nanmean(theta_NREM), nanmean(theta_AW), nanmean(theta_QW)];
            LFPBehavior_CombinedData.AlphaPeriod_dark_baseline = [nanmean(alpha_REM), nanmean(alpha_NREM), nanmean(alpha_AW), nanmean(alpha_QW)];
            LFPBehavior_CombinedData.totalsPeriod_dark_baseline = [total_delta, total_theta, total_alpha];
            LFPBehavior_CombinedData.PercentsPeriod_dark_baseline = [REMtime/total_time, NREMtime/total_time, AWtime/total_time, QWtime/total_time];
            LFPBehavior_CombinedData.DurationsPeriod_dark_baseline = [nanmean(allREMtimes), nanmean(allNREMtimes), nanmean(allAWtimes), nanmean(allQWtimes)];
            
            LFPBehavior_CombinedData.DeltaAnimal_dark_baseline = [];
            LFPBehavior_CombinedData.ThetaAnimal_dark_baseline = [];
            LFPBehavior_CombinedData.AlphaAnimal_dark_baseline = [];
            LFPBehavior_CombinedData.totalsAnimal_dark_baseline = [];
            LFPBehavior_CombinedData.PercentsAnimal_dark_baseline = [];
            LFPBehavior_CombinedData.DurationsAnimal_dark_baseline = [];
            LFPBehavior_CombinedData.ANIMALS_BLdark = {};
        else
            LFPBehavior_CombinedData.DeltaPeriod_dark_baseline(end+1, :) = [nanmean(delta_REM), nanmean(delta_NREM), nanmean(delta_AW), nanmean(delta_QW)];
            LFPBehavior_CombinedData.ThetaPeriod_dark_baseline(end+1, :) = [nanmean(theta_REM), nanmean(theta_NREM), nanmean(theta_AW), nanmean(theta_QW)];
            LFPBehavior_CombinedData.AlphaPeriod_dark_baseline(end+1, :) = [nanmean(alpha_REM), nanmean(alpha_NREM), nanmean(alpha_AW), nanmean(alpha_QW)];
            LFPBehavior_CombinedData.totalsPeriod_dark_baseline(end+1, :) = [total_delta, total_theta, total_alpha];
            LFPBehavior_CombinedData.PercentsPeriod_dark_baseline(end+1, :) = [REMtime/total_time, NREMtime/total_time, AWtime/total_time, QWtime/total_time];
            LFPBehavior_CombinedData.DurationsPeriod_dark_baseline(end+1,:) = [nanmean(allREMtimes), nanmean(allNREMtimes), nanmean(allAWtimes), nanmean(allQWtimes)];
        end
    
    end 
    end
    
    % in case only get one CNO period for an animal, still want to be able
    % to average at least 2 baseline periods because of behavioral
    % variability, but need to record the same number of CNO and baseline
    % periods for stats (paired ttest) purposes, so adjust here
    if CNO_light_num < BL_light_num
        diff = BL_light_num-CNO_light_num;
        LFPBehavior_CombinedData.DeltaPeriod_light_baseline(end-diff, :) = nanmean(LFPBehavior_CombinedData.DeltaPeriod_light_baseline(end-diff:end, :), 1);
        LFPBehavior_CombinedData.DeltaPeriod_light_baseline(end-diff+1:end, :) = [];
        LFPBehavior_CombinedData.ThetaPeriod_light_baseline(end-diff, :) = nanmean(LFPBehavior_CombinedData.ThetaPeriod_light_baseline(end-diff:end, :), 1);
        LFPBehavior_CombinedData.ThetaPeriod_light_baseline(end-diff+1:end, :) = [];
        LFPBehavior_CombinedData.AlphaPeriod_light_baseline(end-diff, :) = nanmean(LFPBehavior_CombinedData.AlphaPeriod_light_baseline(end-diff:end, :), 1);
        LFPBehavior_CombinedData.AlphaPeriod_light_baseline(end-diff+1:end, :) = [];
        LFPBehavior_CombinedData.totalsPeriod_light_baseline(end-diff, :) = nanmean(LFPBehavior_CombinedData.totalsPeriod_light_baseline(end-diff:end, :), 1);
        LFPBehavior_CombinedData.totalsPeriod_light_baseline(end-diff+1:end, :) = [];
        LFPBehavior_CombinedData.PercentsPeriod_light_baseline(end-diff, :) = nanmean(LFPBehavior_CombinedData.PercentsPeriod_light_baseline(end-diff:end, :), 1);
        LFPBehavior_CombinedData.PercentsPeriod_light_baseline(end-diff+1:end, :) = [];
        LFPBehavior_CombinedData.DurationsPeriod_light_baseline(end-diff, :) = nanmean(LFPBehavior_CombinedData.DurationsPeriod_light_baseline(end-diff:end,:), 1);
        LFPBehavior_CombinedData.DurationsPeriod_light_baseline(end-diff+1:end, :) = [];
    end
    
    if CNO_dark_num < BL_dark_num
        diff = BL_dark_num-CNO_dark_num;
        LFPBehavior_CombinedData.DeltaPeriod_dark_baseline(end-diff, :) = nanmean(LFPBehavior_CombinedData.DeltaPeriod_dark_baseline(end-diff:end, :), 1);
        LFPBehavior_CombinedData.DeltaPeriod_dark_baseline(end-diff+1:end, :) = [];
        LFPBehavior_CombinedData.ThetaPeriod_dark_baseline(end-diff, :) = nanmean(LFPBehavior_CombinedData.ThetaPeriod_dark_baseline(end-diff:end, :), 1);
        LFPBehavior_CombinedData.ThetaPeriod_dark_baseline(end-diff+1:end, :) = [];
        LFPBehavior_CombinedData.AlphaPeriod_dark_baseline(end-diff, :) = nanmean(LFPBehavior_CombinedData.AlphaPeriod_dark_baseline(end-diff:end, :), 1);
        LFPBehavior_CombinedData.AlphaPeriod_dark_baseline(end-diff+1:end, :) = [];
        LFPBehavior_CombinedData.totalsPeriod_dark_baseline(end-diff, :) = nanmean(LFPBehavior_CombinedData.totalsPeriod_dark_baseline(end-diff:end, :), 1);
        LFPBehavior_CombinedData.totalsPeriod_dark_baseline(end-diff+1:end, :) = [];
        LFPBehavior_CombinedData.PercentsPeriod_dark_baseline(end-diff, :) = nanmean(LFPBehavior_CombinedData.PercentsPeriod_dark_baseline(end-diff:end, :), 1);
        LFPBehavior_CombinedData.PercentsPeriod_dark_baseline(end-diff+1:end, :) = [];
        LFPBehavior_CombinedData.DurationsPeriod_dark_baseline(end-diff, :) = nanmean(LFPBehavior_CombinedData.DurationsPeriod_dark_baseline(end-diff:end,:), 1);
        LFPBehavior_CombinedData.DurationsPeriod_dark_baseline(end-diff+1:end, :) = [];
    end
    % Now all number of periods should be equal to num of CNO light periods
    % and CNO dark periods
    
    % Store Behavior and LFP info by animal
    
        LFPBehavior_CombinedData.DeltaAnimal_light_CNO(end+1,:) = [nanmean(LFPBehavior_CombinedData.DeltaPeriod_light_CNO(end-(CNO_light_num-1):end,:), 1)];
        LFPBehavior_CombinedData.ThetaAnimal_light_CNO(end+1,:) = [nanmean(LFPBehavior_CombinedData.ThetaPeriod_light_CNO(end-(CNO_light_num-1):end,:), 1)];
        LFPBehavior_CombinedData.AlphaAnimal_light_CNO(end+1,:) = [nanmean(LFPBehavior_CombinedData.AlphaPeriod_light_CNO(end-(CNO_light_num-1):end,:), 1)];
        LFPBehavior_CombinedData.totalsAnimal_light_CNO(end+1,:) = [nanmean(LFPBehavior_CombinedData.totalsPeriod_light_CNO(end-(CNO_light_num-1):end,:), 1)];
        LFPBehavior_CombinedData.PercentsAnimal_light_CNO(end+1,:) = [nanmean(LFPBehavior_CombinedData.PercentsPeriod_light_CNO(end-(CNO_light_num-1):end,:), 1)];
        LFPBehavior_CombinedData.DurationsAnimal_light_CNO(end+1,:) = [nanmean(LFPBehavior_CombinedData.DurationsPeriod_light_CNO(end-(CNO_light_num-1):end,:), 1)];
        LFPBehavior_CombinedData.ANIMALS_CNOlight(end+1:end+CNO_light_num, 1) = {thisAnimal};
    
        LFPBehavior_CombinedData.DeltaAnimal_dark_CNO(end+1,:) = [nanmean(LFPBehavior_CombinedData.DeltaPeriod_dark_CNO(end-(CNO_dark_num-1):end,:), 1)];
        LFPBehavior_CombinedData.ThetaAnimal_dark_CNO(end+1,:) = [nanmean(LFPBehavior_CombinedData.ThetaPeriod_dark_CNO(end-(CNO_dark_num-1):end,:), 1)];
        LFPBehavior_CombinedData.AlphaAnimal_dark_CNO(end+1,:) = [nanmean(LFPBehavior_CombinedData.AlphaPeriod_dark_CNO(end-(CNO_dark_num-1):end,:), 1)];
        LFPBehavior_CombinedData.totalsAnimal_dark_CNO(end+1,:) = [nanmean(LFPBehavior_CombinedData.totalsPeriod_dark_CNO(end-(CNO_dark_num-1):end,:), 1)];
        LFPBehavior_CombinedData.PercentsAnimal_dark_CNO(end+1,:) = [nanmean(LFPBehavior_CombinedData.PercentsPeriod_dark_CNO(end-(CNO_dark_num-1):end,:),1)];
        LFPBehavior_CombinedData.DurationsAnimal_dark_CNO(end+1,:) = [nanmean(LFPBehavior_CombinedData.DurationsPeriod_dark_CNO(end-(CNO_dark_num-1):end,:),1)];
        LFPBehavior_CombinedData.ANIMALS_CNOdark(end+1:end+CNO_dark_num, 1) = {thisAnimal};
    
        LFPBehavior_CombinedData.DeltaAnimal_light_baseline(end+1,:) = [nanmean(LFPBehavior_CombinedData.DeltaPeriod_light_baseline(end-(CNO_light_num-1):end,:), 1)];
        LFPBehavior_CombinedData.ThetaAnimal_light_baseline(end+1,:) = [nanmean(LFPBehavior_CombinedData.ThetaPeriod_light_baseline(end-(CNO_light_num-1):end,:), 1)];
        LFPBehavior_CombinedData.AlphaAnimal_light_baseline(end+1,:) = [nanmean(LFPBehavior_CombinedData.AlphaPeriod_light_baseline(end-(CNO_light_num-1):end,:), 1)];
        LFPBehavior_CombinedData.totalsAnimal_light_baseline(end+1,:) = [nanmean(LFPBehavior_CombinedData.totalsPeriod_light_baseline(end-(CNO_light_num-1):end,:), 1)];
        LFPBehavior_CombinedData.PercentsAnimal_light_baseline(end+1,:) = [nanmean(LFPBehavior_CombinedData.PercentsPeriod_light_baseline(end-(CNO_light_num-1):end,:),1)];
        LFPBehavior_CombinedData.DurationsAnimal_light_baseline(end+1,:) = [nanmean(LFPBehavior_CombinedData.DurationsPeriod_light_baseline(end-(CNO_light_num-1):end,:),1)];
        LFPBehavior_CombinedData.ANIMALS_BLlight(end+1:end+CNO_light_num, 1) = {thisAnimal};
        
        LFPBehavior_CombinedData.DeltaAnimal_dark_baseline(end+1,:) = [nanmean(LFPBehavior_CombinedData.DeltaPeriod_dark_baseline(end-(CNO_dark_num-1):end,:), 1)];
        LFPBehavior_CombinedData.ThetaAnimal_dark_baseline(end+1,:) = [nanmean(LFPBehavior_CombinedData.ThetaPeriod_dark_baseline(end-(CNO_dark_num-1):end,:), 1)];
        LFPBehavior_CombinedData.AlphaAnimal_dark_baseline(end+1,:) = [nanmean(LFPBehavior_CombinedData.AlphaPeriod_dark_baseline(end-(CNO_dark_num-1):end,:), 1)];
        LFPBehavior_CombinedData.totalsAnimal_dark_baseline(end+1,:) = [nanmean(LFPBehavior_CombinedData.totalsPeriod_dark_baseline(end-(CNO_dark_num-1):end,:), 1)];
        LFPBehavior_CombinedData.PercentsAnimal_dark_baseline(end+1,:) = [nanmean(LFPBehavior_CombinedData.PercentsPeriod_dark_baseline(end-(CNO_dark_num-1):end,:),1)];
        LFPBehavior_CombinedData.DurationsAnimal_dark_baseline(end+1,:) = [nanmean(LFPBehavior_CombinedData.DurationsPeriod_dark_baseline(end-(CNO_dark_num-1):end,:),1)];
        LFPBehavior_CombinedData.ANIMALS_BLdark(end+1:end+CNO_dark_num, 1) = {thisAnimal};

end

if pplot == 1
    stats = Plot_LFPBehavior_CNO_LightDarkCombine_JBmanuscript(LFPBehavior_CombinedData);
else
    stats = 0;
end

end