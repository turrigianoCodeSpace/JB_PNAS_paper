%
%This is a script to go through analyzed mini data and compute more
%accurate rise times
% Script was modified from Brian Cary's script by Juliet Bottorff
% ACh manuscript supplemental information (rise/decay times for minis)


% THESE CELLS ALL STARTED DYING/BLOWING UP WITHIN FIRST 20 SWEEPS, SO
% EXCLUDING:
exclude = {'date072221_Cell09',  'date072221_Cell10', 'date072321_Cell04', 'date072321_Cell09', 'date092021_Cell01', 'date092221_Cell09', 'date092221_Cell10', 'date092221_Cell14', 'date092321_Cell04'};

% From loading ALL_MINI_DATA.mat
full_data_strct = analysis_output;

cond1 = 'Hm4di';
cond2 = 'CNOonly';
cond3 = 'BF_Hm4di';
cond4 = 'BF_CNOonly';

cond1_strct = full_data_strct.(char(cond1));
cond2_strct = full_data_strct.(char(cond2));
cond3_strct = full_data_strct.(char(cond3));
cond4_strct = full_data_strct.(char(cond4));

cond1_simple_strct = simple_output.(char(cond1));
cond2_simple_strct = simple_output.(char(cond2));
cond3_simple_strct = simple_output.(char(cond3));
cond4_simple_strct = simple_output.(char(cond4));

% RUN FIRST CONDITION
data_strct = cond1_strct;
cell_names = fields(data_strct);
sample_rate = 10000;

rise_time_col = 14;
shift_axis = 40;

Hm4di_full_taus = [];
Hm4di_full_cell_avgs = {};
Hm4di_full_rise_times = {};
Hm4di_full_cell_rise_times = [];
count = 0;
Hm4di_exclude = [];
for cell_num = 1:length(cell_names)

    if ismember(cond1_simple_strct{cell_num,1}, exclude)
        continue
    end
    
     % exclude cells based on extreme PP numbers
    %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
        thisCell_PP_data = [];
        for col = 2:12
        
        this_entry = cond1_simple_strct{cell_num+1,col};
        this_entry_vect = [];
        
        for entry_cell = 1:length(this_entry)
            if ~isempty(this_entry{entry_cell})
            this_entry_vect(entry_cell) = nanmean(this_entry{entry_cell});
            end
        end
        thisCell_PP_data(col-1) = nanmean(this_entry_vect);
        end
    
        if thisCell_PP_data(2) > 20000000
            continue
        elseif thisCell_PP_data(3) < 90000000
            continue
        elseif thisCell_PP_data(4) < 0.000000000050
            continue
        end
        
        % exclude cells that are too noisy
        cell_RMSE = [];
        for t = 1:size(cond1_strct.(char(cell_names(cell_num))).mini_data, 1)
        TRACE = cond1_strct.(char(cell_names(cell_num))).mini_data{t,2}.DATA;
        av = mean(TRACE);

        RMSE = sqrt(mean((TRACE-av).^2));
        cell_RMSE(end+1) = RMSE;
        end
        
        if mean(cell_RMSE) > 0.0000000000038
            continue
        end
        
        
        
        count = count + 1;
    
    cell_taus = [];
    cell_shift_events = [];
    trace_names = (data_strct.(char(cell_names(cell_num))).mini_data(:,1));
    cell_rise_times = {};
    trace_avg_rise_time = [];
    for trace = 1:length(trace_names)
        trace_shifted_events = [];
        disp(['cell=',num2str(cell_num),' trace=',num2str(trace)])
        timeindx = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.timeindx;
        trace_rise_times = [];
        fprintf('Mini =    ');
        for mini = 1:length(timeindx)
            try
            for char_num = 1:length(num2str(mini))
            fprintf('\b')
            end    
            fprintf('%d',mini);

            start_ind = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.event_start_ind(mini);
            pk_ind = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.sm_pk_ind(mini);
            pk_val = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.amp(mini)*10^-12;
            mini_raw = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.events(:,mini);  

            %%%%%%%% Shift events for alignment
            smooth_mini = smooth(mini_raw(~isnan(mini_raw)),19,'sgolay',3)';

            baseline_cur = nanmedian(smooth_mini(1:start_ind));
            smooth_mini = smooth_mini - baseline_cur;

            end_baseline = nanmean(mini_raw(50:end));

            sm_mini_raw = smooth(mini_raw(~isnan(mini_raw)),19,'sgolay',3)'; % was 13
            sm_mini_raw = sm_mini_raw - baseline_cur;

            near_pk_ind = find(sm_mini_raw < (sm_mini_raw(start_ind)-pk_val)*.80,1);
            shift_thresh = (sm_mini_raw(start_ind)-pk_val)/2;
            [~, half_way_slope] = min(abs(smooth_mini(start_ind:near_pk_ind)-shift_thresh));
            half_way_slope = half_way_slope + start_ind -1;
            shift_ind_2 = shift_axis - half_way_slope -1;

            baseline_mini_raw = (mini_raw(~isnan(mini_raw)) - baseline_cur)';
            if shift_ind_2 > 0 
                shifted_event = [NaN(1,shift_ind_2) baseline_mini_raw];
            elseif shift_ind_2 < 0 
                shifted_event = baseline_mini_raw(-shift_ind_2:end);
            else
                shifted_event = baseline_mini_raw(shift_ind_2+1:end);
            end

            sm_factor = 19; 
            sm_shift_event = smooth(shifted_event,sm_factor,'sgolay',3);  
            if end_baseline>-2e-12
                trace_shifted_events = padmat(trace_shifted_events, sm_shift_event', 1);
            end            
            %%%%%%%%%
            
            smooth_mini = smooth(mini_raw,19); 
            mini_ind = 1:length(smooth_mini);
            inter_time = 1:0.1:length(smooth_mini);
            interp_sm_mini = interp1(mini_ind,smooth_mini,inter_time);

            baseline_cur = nanmedian(smooth_mini(1:start_ind));

            pk_val = smooth_mini(pk_ind);
            amp_rise = abs(pk_val) - abs(baseline_cur);
            
            flp_mini = flip(interp_sm_mini(1:pk_ind*10));
            
            rise_t1_flp = find(flp_mini>=(baseline_cur-amp_rise*0.1),1);
            rise_t2_flp = find(flp_mini>=(-abs(pk_val)+amp_rise*0.1),1);
            
            rise_t1 = pk_ind*10-rise_t1_flp;
            rise_t2 = pk_ind*10-rise_t2_flp;
            
            if isempty(rise_t2)
                rise_t2 = pk_ind*10;
            end
            trace_rise_times(mini) = (rise_t2-rise_t1)/(sample_rate*10);

            catch
                trace_rise_times(mini) = NaN;
                disp('err')
            end

            
        end
        fprintf('\n');
        cell_rise_times{trace,1} = trace_rise_times;
        trace_avg_rise_time(end+1) = nanmean(trace_rise_times);
        cell_shift_events = padmat(cell_shift_events,trace_shifted_events,1);
    end
    cell_avg_rise_time = nanmean(trace_avg_rise_time);
    
    cell_mean_trace = nanmean(cell_shift_events',2);
    
    pk_indx_st = 20;
    peak_wind_indx = pk_indx_st:100;
    length_fit_indx = 150;
    
    [min_val, min_ind] = min(cell_mean_trace(peak_wind_indx));
    pk_ind = min_ind+pk_indx_st;
    
    try
        fit_data = cell_mean_trace(pk_ind:pk_ind+length_fit_indx);
    catch
        fit_data = cell_mean_trace(pk_ind:end);
    end
    
    fit_data = fit_data.*10^12;
    tau_est = 0.003;
    time = (0:(1/sample_rate):(length(fit_data)/sample_rate)-(1/sample_rate)).';
    s = fitoptions('Method', 'NonlinearLeastSquares', ...
        'StartPoint', [mean(fit_data(1:5)), tau_est],...
         'Lower', [mean(fit_data(1:3))*1.1, 0.001],...
         'Upper', [mean(fit_data(1:3))*0.9, 0.006]);
    f = fittype('a*(exp(-x/b))','options',s);

    [exp_fit,gof] = fit(time,fit_data,f);
    cval = coeffvalues(exp_fit);

    Hm4di_full_taus(count) = cval(2);
    if cval(2) > 0.004
        Hm4di_exclude = [Hm4di_exclude, count];
    end
    Hm4di_full_cell_avgs{count,1} = cell_mean_trace;
    Hm4di_full_rise_times{count,1} = cell_rise_times;
    Hm4di_full_cell_rise_times(count) = cell_avg_rise_time;
end


% RUN SECOND CONDITION
data_strct = cond2_strct;
cell_names = fields(data_strct);
sample_rate = 10000;

rise_time_col = 14;
shift_axis = 40;

CNOonly_full_taus = [];
CNOonly_full_cell_avgs = {};
%cell_shift_events = [];
CNOonly_full_rise_times = {};
CNOonly_full_cell_rise_times = [];
count = 0;
CNOonly_exclude = [];
for cell_num = 1:length(cell_names)
    
    if ismember(cond2_simple_strct{cell_num,1}, exclude)
        continue
    end
    
    
    % exclude cells based on extreme PP numbers
    %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
        thisCell_PP_data = [];
        for col = 2:12
        
        this_entry = cond2_simple_strct{cell_num+1,col};
        this_entry_vect = [];
        
        for entry_cell = 1:length(this_entry)
            if ~isempty(this_entry{entry_cell})
            this_entry_vect(entry_cell) = nanmean(this_entry{entry_cell});
            end
        end
        thisCell_PP_data(col-1) = nanmean(this_entry_vect);
        end
    
        if thisCell_PP_data(2) > 20000000
            continue
        elseif thisCell_PP_data(3) < 90000000 
            continue
        elseif thisCell_PP_data(4) < 0.000000000050
            continue
        end
        
        % exclude cells that are too noisy
        cell_RMSE = [];
        for t = 1:size(cond2_strct.(char(cell_names(cell_num))).mini_data, 1)
        TRACE = cond2_strct.(char(cell_names(cell_num))).mini_data{t,2}.DATA;
        av = mean(TRACE);

        RMSE = sqrt(mean((TRACE-av).^2));
        cell_RMSE(end+1) = RMSE;
        end
        
        if mean(cell_RMSE) > 0.0000000000038
            continue
        end
        
        count = count + 1;
        
    
    cell_taus = [];
    cell_shift_events = [];
    trace_names = (data_strct.(char(cell_names(cell_num))).mini_data(:,1));
    cell_rise_times = {};
    trace_avg_rise_time = [];
    for trace = 1:length(trace_names)
        trace_shifted_events = [];
        disp(['cell=',num2str(cell_num),' trace=',num2str(trace)])
        timeindx = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.timeindx;
        trace_rise_times = [];
        fprintf('Mini =    ');
        for mini = 1:length(timeindx)
            try
            for char_num = 1:length(num2str(mini))
            fprintf('\b')
            end    
            fprintf('%d',mini);

            start_ind = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.event_start_ind(mini);
            pk_ind = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.sm_pk_ind(mini);
            pk_val = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.amp(mini)*10^-12;
            mini_raw = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.events(:,mini);  

            %%%%%%%% Shift events for alignment
            smooth_mini = smooth(mini_raw(~isnan(mini_raw)),19,'sgolay',3)';

            baseline_cur = nanmedian(smooth_mini(1:start_ind));
            smooth_mini = smooth_mini - baseline_cur;

            end_baseline = nanmean(mini_raw(50:end));

            sm_mini_raw = smooth(mini_raw(~isnan(mini_raw)),19,'sgolay',3)'; % was 13
            sm_mini_raw = sm_mini_raw - baseline_cur;

            near_pk_ind = find(sm_mini_raw < (sm_mini_raw(start_ind)-pk_val)*.80,1);
            shift_thresh = (sm_mini_raw(start_ind)-pk_val)/2;
            [~, half_way_slope] = min(abs(smooth_mini(start_ind:near_pk_ind)-shift_thresh));
            half_way_slope = half_way_slope + start_ind -1;
            shift_ind_2 = shift_axis - half_way_slope -1;

            baseline_mini_raw = (mini_raw(~isnan(mini_raw)) - baseline_cur)';
            if shift_ind_2 > 0 
                shifted_event = [NaN(1,shift_ind_2) baseline_mini_raw];
            elseif shift_ind_2 < 0 
                shifted_event = baseline_mini_raw(-shift_ind_2:end);
            else
                shifted_event = baseline_mini_raw(shift_ind_2+1:end);
            end

            sm_factor = 19;
            sm_shift_event = smooth(shifted_event,sm_factor,'sgolay',3); 
            if end_baseline>-2e-12
                trace_shifted_events = padmat(trace_shifted_events, sm_shift_event', 1);
            end            
            %%%%%%%%%
            
            smooth_mini = smooth(mini_raw,19); 
            mini_ind = 1:length(smooth_mini);
            inter_time = 1:0.1:length(smooth_mini);
            interp_sm_mini = interp1(mini_ind,smooth_mini,inter_time);

            baseline_cur = nanmedian(smooth_mini(1:start_ind));

            pk_val = smooth_mini(pk_ind);
            amp_rise = abs(pk_val) - abs(baseline_cur);
            
            flp_mini = flip(interp_sm_mini(1:pk_ind*10));
            
            rise_t1_flp = find(flp_mini>=(baseline_cur-amp_rise*0.1),1);
            rise_t2_flp = find(flp_mini>=(-abs(pk_val)+amp_rise*0.1),1);
            
            rise_t1 = pk_ind*10-rise_t1_flp;
            rise_t2 = pk_ind*10-rise_t2_flp;
            
            if isempty(rise_t2)
                rise_t2 = pk_ind*10;
            end
            trace_rise_times(mini) = (rise_t2-rise_t1)/(sample_rate*10);
            

            catch
                trace_rise_times(mini) = NaN;
                disp('err')
            end
            
        end
        fprintf('\n');
        cell_rise_times{trace,1} = trace_rise_times;
        trace_avg_rise_time(end+1) = nanmean(trace_rise_times);
        cell_shift_events = padmat(cell_shift_events,trace_shifted_events,1);
    end
    cell_avg_rise_time = nanmean(trace_avg_rise_time);
    
    cell_mean_trace = nanmean(cell_shift_events',2);

    pk_indx_st = 20;
    peak_wind_indx = pk_indx_st:100;
    length_fit_indx = 150;
    
    [min_val, min_ind] = min(cell_mean_trace(peak_wind_indx));
    pk_ind = min_ind+pk_indx_st;
    
    try
        fit_data = cell_mean_trace(pk_ind:pk_ind+length_fit_indx);
    catch
        fit_data = cell_mean_trace(pk_ind:end);
    end
    
    fit_data = fit_data.*10^12;
    tau_est = 0.003;
    time = (0:(1/sample_rate):(length(fit_data)/sample_rate)-(1/sample_rate)).';
    s = fitoptions('Method', 'NonlinearLeastSquares', ...
        'StartPoint', [mean(fit_data(1:5)), tau_est],...
         'Lower', [mean(fit_data(1:3))*1.1, 0.001],...
         'Upper', [mean(fit_data(1:3))*0.9, 0.006]);
    f = fittype('a*(exp(-x/b))','options',s);

    [exp_fit,gof] = fit(time,fit_data,f);
    cval = coeffvalues(exp_fit);

    CNOonly_full_taus(count) = cval(2);
    if cval(2) > 0.004
        CNOonly_exclude = [CNOonly_exclude, count];
    end
    CNOonly_full_cell_avgs{count,1} = cell_mean_trace;
    CNOonly_full_rise_times{count,1} = cell_rise_times;
    CNOonly_full_cell_rise_times(count) = cell_avg_rise_time;
end

% RUN THIRD CONDITION
data_strct = cond3_strct;
cell_names = fields(data_strct);
sample_rate = 10000;

rise_time_col = 14;
shift_axis = 40;

BF_Hm4di_full_taus = [];
BF_Hm4di_full_cell_avgs = {};
%cell_shift_events = [];
BF_Hm4di_full_rise_times = {};
BF_Hm4di_full_cell_rise_times = [];
count = 0;
BF_Hm4di_exclude = [];
for cell_num = 1:length(cell_names)
    
    if ismember(cond3_simple_strct{cell_num,1}, exclude)
        continue
    end
    
    
    % exclude cells based on extreme PP numbers
    %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
        thisCell_PP_data = [];
        for col = 2:12
        
        this_entry = cond3_simple_strct{cell_num+1,col};
        this_entry_vect = [];
        
        for entry_cell = 1:length(this_entry)
            if ~isempty(this_entry{entry_cell})
            this_entry_vect(entry_cell) = nanmean(this_entry{entry_cell});
            end
        end
        thisCell_PP_data(col-1) = nanmean(this_entry_vect);
        end
    
        if thisCell_PP_data(2) > 20000000
            continue
        elseif thisCell_PP_data(3) < 90000000 
            continue
        elseif thisCell_PP_data(4) < 0.000000000050
            continue
        end
        
        % exclude cells that are too noisy
        cell_RMSE = [];
        for t = 1:size(cond3_strct.(char(cell_names(cell_num))).mini_data, 1)
        TRACE = cond3_strct.(char(cell_names(cell_num))).mini_data{t,2}.DATA;
        av = mean(TRACE);

        RMSE = sqrt(mean((TRACE-av).^2));
        cell_RMSE(end+1) = RMSE;
        end
        
        if mean(cell_RMSE) > 0.0000000000038
            continue
        end
        
        count = count + 1;
        
    
    cell_taus = [];
    cell_shift_events = [];
    trace_names = (data_strct.(char(cell_names(cell_num))).mini_data(:,1));
    cell_rise_times = {};
    trace_avg_rise_time = [];
    for trace = 1:length(trace_names)
        trace_shifted_events = [];
        disp(['cell=',num2str(cell_num),' trace=',num2str(trace)])
        timeindx = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.timeindx;
        trace_rise_times = [];
        fprintf('Mini =    ');
        for mini = 1:length(timeindx)
            try
            for char_num = 1:length(num2str(mini))
            fprintf('\b')
            end    
            fprintf('%d',mini);

            start_ind = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.event_start_ind(mini);
            pk_ind = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.sm_pk_ind(mini);
            pk_val = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.amp(mini)*10^-12;
            mini_raw = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.events(:,mini);  

            %%%%%%%% Shift events for alignment
            smooth_mini = smooth(mini_raw(~isnan(mini_raw)),19,'sgolay',3)';

            baseline_cur = nanmedian(smooth_mini(1:start_ind));
            smooth_mini = smooth_mini - baseline_cur;

            end_baseline = nanmean(mini_raw(50:end));

            sm_mini_raw = smooth(mini_raw(~isnan(mini_raw)),19,'sgolay',3)'; % was 13
            sm_mini_raw = sm_mini_raw - baseline_cur;

            near_pk_ind = find(sm_mini_raw < (sm_mini_raw(start_ind)-pk_val)*.80,1);
            shift_thresh = (sm_mini_raw(start_ind)-pk_val)/2;
            [~, half_way_slope] = min(abs(smooth_mini(start_ind:near_pk_ind)-shift_thresh));
            half_way_slope = half_way_slope + start_ind -1;
            shift_ind_2 = shift_axis - half_way_slope -1;

            baseline_mini_raw = (mini_raw(~isnan(mini_raw)) - baseline_cur)';
            if shift_ind_2 > 0 
                shifted_event = [NaN(1,shift_ind_2) baseline_mini_raw];
            elseif shift_ind_2 < 0 
                shifted_event = baseline_mini_raw(-shift_ind_2:end);
            else
                shifted_event = baseline_mini_raw(shift_ind_2+1:end);
            end

            sm_factor = 19; 
            sm_shift_event = smooth(shifted_event,sm_factor,'sgolay',3); 
            if end_baseline>-2e-12
                trace_shifted_events = padmat(trace_shifted_events, sm_shift_event', 1);
            end            
            %%%%%%%%%
            
            smooth_mini = smooth(mini_raw,19); 
            mini_ind = 1:length(smooth_mini);
            inter_time = 1:0.1:length(smooth_mini);
            interp_sm_mini = interp1(mini_ind,smooth_mini,inter_time);

            baseline_cur = nanmedian(smooth_mini(1:start_ind));

            pk_val = smooth_mini(pk_ind);
            amp_rise = abs(pk_val) - abs(baseline_cur);
            
            flp_mini = flip(interp_sm_mini(1:pk_ind*10));
            
            rise_t1_flp = find(flp_mini>=(baseline_cur-amp_rise*0.1),1);
            rise_t2_flp = find(flp_mini>=(-abs(pk_val)+amp_rise*0.1),1);
            
            rise_t1 = pk_ind*10-rise_t1_flp;
            rise_t2 = pk_ind*10-rise_t2_flp;
            
            if isempty(rise_t2)
                rise_t2 = pk_ind*10;
            end
            trace_rise_times(mini) = (rise_t2-rise_t1)/(sample_rate*10);
            

            catch
                trace_rise_times(mini) = NaN;
                disp('err')
            end

            
        end
        fprintf('\n');
        cell_rise_times{trace,1} = trace_rise_times;
        trace_avg_rise_time(end+1) = nanmean(trace_rise_times);
        cell_shift_events = padmat(cell_shift_events,trace_shifted_events,1);
    end
    cell_avg_rise_time = nanmean(trace_avg_rise_time);
    
    cell_mean_trace = nanmean(cell_shift_events',2);
lim([0, 150])
    
    pk_indx_st = 20;
    peak_wind_indx = pk_indx_st:100;
    length_fit_indx = 150;
    
    [min_val, min_ind] = min(cell_mean_trace(peak_wind_indx));
    pk_ind = min_ind+pk_indx_st;
    
    try
        fit_data = cell_mean_trace(pk_ind:pk_ind+length_fit_indx);
    catch
        fit_data = cell_mean_trace(pk_ind:end);
    end
    
    fit_data = fit_data.*10^12;
    tau_est = 0.003;
    time = (0:(1/sample_rate):(length(fit_data)/sample_rate)-(1/sample_rate)).';
    s = fitoptions('Method', 'NonlinearLeastSquares', ...
        'StartPoint', [mean(fit_data(1:5)), tau_est],...
         'Lower', [mean(fit_data(1:3))*1.1, 0.001],...
         'Upper', [mean(fit_data(1:3))*0.9, 0.006]);
    f = fittype('a*(exp(-x/b))','options',s);

    [exp_fit,gof] = fit(time,fit_data,f);
    cval = coeffvalues(exp_fit);

    BF_Hm4di_full_taus(count) = cval(2);
    if cval(2) > 0.004
        BF_Hm4di_exclude = [BF_Hm4di_exclude, count];
    end
    BF_Hm4di_full_cell_avgs{count,1} = cell_mean_trace;
    BF_Hm4di_full_rise_times{count,1} = cell_rise_times;
    BF_Hm4di_full_cell_rise_times(count) = cell_avg_rise_time;
end

% RUN FOURTH CONDITION
data_strct = cond4_strct;
cell_names = fields(data_strct);
sample_rate = 10000;

rise_time_col = 14;
shift_axis = 40;

BF_CNOonly_full_taus = [];
BF_CNOonly_full_cell_avgs = {};
BF_CNOonly_full_rise_times = {};
BF_CNOonly_full_cell_rise_times = [];
count = 0;
BF_CNOonly_exclude = [];
for cell_num = 1:length(cell_names)
    
    if ismember(cond4_simple_strct{cell_num,1}, exclude)
        continue
    end
    
    
    % exclude cells based on extreme PP numbers
    %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
        thisCell_PP_data = [];
        for col = 2:12
        
        this_entry = cond4_simple_strct{cell_num+1,col};
        this_entry_vect = [];
        
        for entry_cell = 1:length(this_entry)
            if ~isempty(this_entry{entry_cell})
            this_entry_vect(entry_cell) = nanmean(this_entry{entry_cell});
            end
        end
        thisCell_PP_data(col-1) = nanmean(this_entry_vect);
        end
    
        if thisCell_PP_data(2) > 20000000
            continue
        elseif thisCell_PP_data(3) < 90000000 
            continue
        elseif thisCell_PP_data(4) < 0.000000000050
            continue
        end
        
        % exclude cells that are too noisy
        cell_RMSE = [];
        for t = 1:size(cond4_strct.(char(cell_names(cell_num))).mini_data, 1)
        TRACE = cond4_strct.(char(cell_names(cell_num))).mini_data{t,2}.DATA;
        av = mean(TRACE);

        RMSE = sqrt(mean((TRACE-av).^2));
        cell_RMSE(end+1) = RMSE;
        end
        
        if mean(cell_RMSE) > 0.0000000000038
            continue
        end
        
        count = count + 1;
        
    
    cell_taus = [];
    cell_shift_events = [];
    trace_names = (data_strct.(char(cell_names(cell_num))).mini_data(:,1));
    cell_rise_times = {};
    trace_avg_rise_time = [];
    for trace = 1:length(trace_names)
        trace_shifted_events = [];
        disp(['cell=',num2str(cell_num),' trace=',num2str(trace)])
        timeindx = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.timeindx;
        trace_rise_times = [];
        fprintf('Mini =    ');
        for mini = 1:length(timeindx)
            try
            for char_num = 1:length(num2str(mini))
            fprintf('\b')
            end    
            fprintf('%d',mini);

            start_ind = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.event_start_ind(mini);
            pk_ind = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.sm_pk_ind(mini);
            pk_val = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.amp(mini)*10^-12;
            mini_raw = data_strct.(char(cell_names(cell_num))).mini_data{trace, 2}.events(:,mini);  

            %%%%%%%% Shift events for alignment
            smooth_mini = smooth(mini_raw(~isnan(mini_raw)),19,'sgolay',3)';

            baseline_cur = nanmedian(smooth_mini(1:start_ind));
            smooth_mini = smooth_mini - baseline_cur;

            end_baseline = nanmean(mini_raw(50:end));

            sm_mini_raw = smooth(mini_raw(~isnan(mini_raw)),19,'sgolay',3)'; % was 13
            sm_mini_raw = sm_mini_raw - baseline_cur;

            near_pk_ind = find(sm_mini_raw < (sm_mini_raw(start_ind)-pk_val)*.80,1);
            shift_thresh = (sm_mini_raw(start_ind)-pk_val)/2;
            [~, half_way_slope] = min(abs(smooth_mini(start_ind:near_pk_ind)-shift_thresh));
            half_way_slope = half_way_slope + start_ind -1;
            shift_ind_2 = shift_axis - half_way_slope -1;

            baseline_mini_raw = (mini_raw(~isnan(mini_raw)) - baseline_cur)';
            if shift_ind_2 > 0 
                shifted_event = [NaN(1,shift_ind_2) baseline_mini_raw];
            elseif shift_ind_2 < 0 
                shifted_event = baseline_mini_raw(-shift_ind_2:end);
            else
                shifted_event = baseline_mini_raw(shift_ind_2+1:end);
            end

            sm_factor = 19; 
            sm_shift_event = smooth(shifted_event,sm_factor,'sgolay',3); 
            if end_baseline>-2e-12
                trace_shifted_events = padmat(trace_shifted_events, sm_shift_event', 1);
            end            
            %%%%%%%%%
            
            smooth_mini = smooth(mini_raw,19);
            mini_ind = 1:length(smooth_mini);
            inter_time = 1:0.1:length(smooth_mini);
            interp_sm_mini = interp1(mini_ind,smooth_mini,inter_time);

            baseline_cur = nanmedian(smooth_mini(1:start_ind));

            pk_val = smooth_mini(pk_ind);
            amp_rise = abs(pk_val) - abs(baseline_cur);
            
            flp_mini = flip(interp_sm_mini(1:pk_ind*10));
            
            rise_t1_flp = find(flp_mini>=(baseline_cur-amp_rise*0.1),1);
            rise_t2_flp = find(flp_mini>=(-abs(pk_val)+amp_rise*0.1),1);
            
            rise_t1 = pk_ind*10-rise_t1_flp;
            rise_t2 = pk_ind*10-rise_t2_flp;
            
            if isempty(rise_t2)
                rise_t2 = pk_ind*10;
            end
            trace_rise_times(mini) = (rise_t2-rise_t1)/(sample_rate*10);
            

            catch
                trace_rise_times(mini) = NaN;
                disp('err')
            end

            
        end
        fprintf('\n');
        cell_rise_times{trace,1} = trace_rise_times;
        trace_avg_rise_time(end+1) = nanmean(trace_rise_times);
        cell_shift_events = padmat(cell_shift_events,trace_shifted_events,1);
    end
    cell_avg_rise_time = nanmean(trace_avg_rise_time);
    
    cell_mean_trace = nanmean(cell_shift_events',2);
 
    pk_indx_st = 20;
    peak_wind_indx = pk_indx_st:100;
    length_fit_indx = 150;
    
    [min_val, min_ind] = min(cell_mean_trace(peak_wind_indx));
    pk_ind = min_ind+pk_indx_st;
    
    try
        fit_data = cell_mean_trace(pk_ind:pk_ind+length_fit_indx);
    catch
        fit_data = cell_mean_trace(pk_ind:end);
    end
    
    fit_data = fit_data.*10^12;
    tau_est = 0.003;
    time = (0:(1/sample_rate):(length(fit_data)/sample_rate)-(1/sample_rate)).';
    s = fitoptions('Method', 'NonlinearLeastSquares', ...
        'StartPoint', [mean(fit_data(1:5)), tau_est],...
         'Lower', [mean(fit_data(1:3))*1.1, 0.001],...
         'Upper', [mean(fit_data(1:3))*0.9, 0.006]);
    f = fittype('a*(exp(-x/b))','options',s);

    [exp_fit,gof] = fit(time,fit_data,f);
    cval = coeffvalues(exp_fit);

    BF_CNOonly_full_taus(count) = cval(2);
    if cval(2) > 0.004
        BF_CNOonly_exclude = [BF_CNOonly_exclude, count];
    end
    BF_CNOonly_full_cell_avgs{count,1} = cell_mean_trace;
    BF_CNOonly_full_rise_times{count,1} = cell_rise_times;
    BF_CNOonly_full_cell_rise_times(count) = cell_avg_rise_time;
end

Colors = [0 0 0; 1 0 1; 0.5 0.5 0.5; 0.58 0.82 0.98];

%risetimes

figure(100);
deviation = 0.35;
Hm4di_x_jitter = 2 + rand(1,length(Hm4di_full_cell_rise_times)).*deviation - deviation/2;
CNOonly_x_jitter = 1 + rand(1,length(CNOonly_full_cell_rise_times)).*deviation - deviation/2;
BF_Hm4di_x_jitter = 4 + rand(1,length(BF_Hm4di_full_cell_rise_times)).*deviation - deviation/2;
BF_CNOonly_x_jitter = 3 + rand(1,length(BF_CNOonly_full_cell_rise_times)).*deviation - deviation/2;

plot(Hm4di_x_jitter,Hm4di_full_cell_rise_times,'o', 'Color', Colors(2,:))
hold on
plot(CNOonly_x_jitter,CNOonly_full_cell_rise_times,'o', 'Color', Colors(1,:))
plot(BF_CNOonly_x_jitter,BF_CNOonly_full_cell_rise_times,'o', 'Color', Colors(3,:))
plot(BF_Hm4di_x_jitter,BF_Hm4di_full_cell_rise_times,'o', 'Color', Colors(4,:))

xlim([0 5])
ylim([0 2.5e-3])
bar(2,nanmean(Hm4di_full_cell_rise_times),'FaceColor','none')
bar(1,nanmean(CNOonly_full_cell_rise_times),'FaceColor','none')
bar(4,nanmean(BF_Hm4di_full_cell_rise_times),'FaceColor','none')
bar(3,nanmean(BF_CNOonly_full_cell_rise_times),'FaceColor','none')
title('rise time');
set(gca,'XTickLabel',{'CNO only', 'V1 Hm4di + CNO', 'BF ACh Hm4di + CNO', 'BF ACh Hm4di + V1 Hm4di + CNO'}, 'box', 'off', 'FontSize', 15)
xticks([1 2 3 4])
xtickangle(45)

rise_data_1 = padmat(CNOonly_full_cell_rise_times', Hm4di_full_cell_rise_times', 2);
rise_data_2 = padmat(BF_CNOonly_full_cell_rise_times', BF_Hm4di_full_cell_rise_times', 2);
rise_data = padmat(rise_data_1, rise_data_2, 2);
[p_kruskalwallis_rise, tbl, stats] = kruskalwallis(rise_data)
rise_stats = multcompare(stats)


% TAUS

figure(200);
deviation = 0.35;
Hm4di_x_jitter = 2 + rand(1,length(Hm4di_full_taus)).*deviation - deviation/2;
CNOonly_x_jitter = 1 + rand(1,length(CNOonly_full_taus)).*deviation - deviation/2;
BF_Hm4di_x_jitter = 4 + rand(1,length(BF_Hm4di_full_taus)).*deviation - deviation/2;
BF_CNOonly_x_jitter = 3 + rand(1,length(BF_CNOonly_full_taus)).*deviation - deviation/2;

plot(Hm4di_x_jitter,Hm4di_full_taus,'o', 'Color', Colors(2,:))
hold on
plot(CNOonly_x_jitter,CNOonly_full_taus,'o', 'Color', Colors(1,:))
plot(BF_Hm4di_x_jitter,BF_Hm4di_full_taus,'o', 'Color', Colors(4,:))
plot(BF_CNOonly_x_jitter,BF_CNOonly_full_taus,'o', 'Color', Colors(3,:))
xlim([0 5])
ylim([0 6.5e-3])
bar(2,nanmean(Hm4di_full_taus),'FaceColor','none')
bar(1,nanmean(CNOonly_full_taus),'FaceColor','none')
bar(4,nanmean(BF_Hm4di_full_taus),'FaceColor','none')
bar(3,nanmean(BF_CNOonly_full_taus),'FaceColor','none')
title('t decay');
set(gca,'XTickLabel',{'CNO only', 'V1 Hm4di + CNO', 'BF ACh Hm4di + CNO' 'BF ACh Hm4di + V1 Hm4di + CNO'}, 'box', 'off', 'FontSize', 15)
xticks([1 2 3 4])
xtickangle(45)

decay_data_1 = padmat(CNOonly_full_taus', Hm4di_full_taus', 2);
decay_data_2 = padmat(BF_CNOonly_full_taus', BF_Hm4di_full_taus', 2);
decay_data = padmat(decay_data_1, decay_data_2, 2);
[p_kruskalwallis_decay, tbl, stats] = kruskalwallis(decay_data)
decay_stats = multcompare(stats)

%plot avg waveforms

Hm4di_wave_avgs =[];
CNOonly_wave_avgs = [];
BF_Hm4di_wave_avgs =[];
BF_CNOonly_wave_avgs = [];

for cell_num = 1:length(Hm4di_full_cell_avgs)
    Hm4di_wave_avgs = padmat(Hm4di_wave_avgs,Hm4di_full_cell_avgs{cell_num},2);
end
for cell_num = 1:length(CNOonly_full_cell_avgs)
    CNOonly_wave_avgs = padmat(CNOonly_wave_avgs,CNOonly_full_cell_avgs{cell_num},2);
end
for cell_num = 1:length(BF_Hm4di_full_cell_avgs)
    BF_Hm4di_wave_avgs = padmat(BF_Hm4di_wave_avgs,BF_Hm4di_full_cell_avgs{cell_num},2);
end
for cell_num = 1:length(BF_CNOonly_full_cell_avgs)
    BF_CNOonly_wave_avgs = padmat(BF_CNOonly_wave_avgs,BF_CNOonly_full_cell_avgs{cell_num},2);
end

mean_Hm4di = nanmean(Hm4di_wave_avgs,2);
mean_CNOonly = nanmean(CNOonly_wave_avgs,2);
mean_BF_Hm4di = nanmean(BF_Hm4di_wave_avgs,2);
mean_BF_CNOonly = nanmean(BF_CNOonly_wave_avgs,2);
sample_rate = 10000;

time = (0:(1/sample_rate):(length(mean_CNOonly)/sample_rate)-(1/sample_rate)).';
figure
plot(time,CNOonly_wave_avgs,'Color',[0 0 0 0.25]); hold on
plot(time,mean_CNOonly,'k','LineWidth',4.0)
xlim([0.0028 0.015])
ylim([-16e-12 2e-12])
title ('CNO only cells - average traces')


time = (0:(1/sample_rate):(length(mean_Hm4di)/sample_rate)-(1/sample_rate)).';
figure
plot(time,Hm4di_wave_avgs,'Color',[0 0 0 0.25]); hold on
plot(time,mean_Hm4di,'m','LineWidth',4.0)
xlim([0.0028 0.015])
ylim([-16e-12 2e-12])
title ('Hm4di cells - average traces')

time = (0:(1/sample_rate):(length(mean_BF_Hm4di)/sample_rate)-(1/sample_rate)).';
figure
plot(time,BF_Hm4di_wave_avgs,'Color',[0 0 0 0.25]); hold on
plot(time,mean_BF_Hm4di,'b','LineWidth',4.0)
xlim([0.0028 0.015])
ylim([-16e-12 2e-12])
title ('BF Hm4di cells - average traces')

time = (0:(1/sample_rate):(length(mean_BF_CNOonly)/sample_rate)-(1/sample_rate)).';
figure
plot(time,BF_CNOonly_wave_avgs,'Color',[0 0 0 0.25]); hold on
plot(time,mean_BF_CNOonly,'k','LineWidth',4.0)
xlim([0.0028 0.015])
ylim([-16e-12 2e-12])
title ('BF CNOonly cells - average traces')


figure_coords = [412 655 200 331];
figure('Position',figure_coords)
hold on
plot(mean_Hm4di,'Color', Colors(2,:),'LineWidth',4.0)

plot([100-55 100-5], [-12e-12 -12e-12], 'k', 'LineWidth', 4) % 5 ms, 5 pA
plot([100-5 100-5], [-7e-12 -12e-12], 'k', 'LineWidth', 4)

box off
xlim([20 130]) 
ylim([-16e-12 2e-12])

figure_coords = [412 655 200 331];
figure('Position',figure_coords)
hold on
plot(mean_CNOonly,'Color', Colors(1,:),'LineWidth',4.0)

plot([100-55 100-5], [-12e-12 -12e-12], 'k', 'LineWidth', 4) % 5 ms, 5 pA
plot([100-5 100-5], [-7e-12 -12e-12], 'k', 'LineWidth', 4)

box off
xlim([20 130]) 
ylim([-16e-12 2e-12])

figure_coords = [412 655 200 331];
figure('Position',figure_coords)
hold on
plot(mean_BF_Hm4di,'Color', Colors(4,:),'LineWidth',4.0)

plot([100-55 100-5], [-12e-12 -12e-12], 'k', 'LineWidth', 4) % 5 ms, 5 pA
plot([100-5 100-5], [-7e-12 -12e-12], 'k', 'LineWidth', 4)

box off
xlim([20 130]) 
ylim([-16e-12 2e-12])

figure_coords = [412 655 200 331];
figure('Position',figure_coords)
hold on
plot(mean_BF_CNOonly,'Color', Colors(3,:),'LineWidth',4.0)

plot([100-55 100-5], [-12e-12 -12e-12], 'k', 'LineWidth', 4) % 5 ms, 5 pA
plot([100-5 100-5], [-7e-12 -12e-12], 'k', 'LineWidth', 4)

box off
xlim([20 130]) 
ylim([-16e-12 2e-12])


pk_indx_st = 20;
peak_wind_indx = pk_indx_st:100;

CNOonly_sc_avg = -mean_CNOonly./min(mean_CNOonly(peak_wind_indx));
Hm4di_sc_avg = -mean_Hm4di./min(mean_Hm4di(peak_wind_indx));
BF_CNOonly_sc_avg = -mean_BF_CNOonly./min(mean_BF_CNOonly(peak_wind_indx));
BF_Hm4di_sc_avg = -mean_BF_Hm4di./min(mean_BF_Hm4di(peak_wind_indx));

figure
hold on
plot(CNOonly_sc_avg, 'Color', Colors(1,:),'LineWidth',2.0)
plot(Hm4di_sc_avg,'Color', Colors(2,:), 'LineWidth',2.0)
plot(BF_CNOonly_sc_avg, 'Color', Colors(3,:),'LineWidth',2.0)
plot(BF_Hm4di_sc_avg,'Color', Colors(4,:), 'LineWidth',2.0)

xlim([30 150])
legend('CNO only', 'V1 Hm4di + CNO', 'BF ACh Hm4di + CNO', 'BF ACh Hm4di + V1 Hm4di + CNO')


