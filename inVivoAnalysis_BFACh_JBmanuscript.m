% FR plots for BF ACh inhib / CNO only datasets; Figure 6E,F
% Bottorff manuscript - updated Dec. 2022


function CELL_DATA = inVivoAnalysis_BFACh_JBmanuscript(MDCELLS, plot_on)

%% filter to only get quality 1,2 cells
all_quals = [MDCELLS.quality];

quality_1_2=find(all_quals<3); 
MDCELLS = MDCELLS(quality_1_2);

%%

%Pick out cells that are active for some percent of total recording time to
%find continuous cells
percent_ON = .7;

expt_end=7; % how many days did the experiment last

lastcell=length(MDCELLS);
ON_cells=0; % to keep track of how many continuous ('on') cells we have
ON_cells_index=NaN(1,lastcell); % initializing index in original matrix of cells that are found to be continuous

for n=1:lastcell
    
    experiment_time_start=MDCELLS(n).trem; % in seconds, when this cell's experiment started
    experiment_time_end=experiment_time_start+(3600*24*expt_end); % in seconds
    total_experiment_time=experiment_time_end-experiment_time_start; % for this cell specifically, total possible time to be continuous or 'ON'
    
        onTimes=MDCELLS(n).onTime;
        offTimes=MDCELLS(n).offTime;
        % don't evaluate cells with mismatched on/off times - need to fix!
        if length(onTimes)~=length(offTimes)
            disp ('This cell has mismatched on/off times - skipping')
            disp (n)
            continue
        end
        
        % loop through all on/off times to find total length of time in seconds that cell is 'on'
        all_ontimes=offTimes(1)-onTimes(1); 
        if length(onTimes)>1 

            for i=2:length(onTimes)
                all_ontimes=all_ontimes+(offTimes(i)-onTimes(i)); 
            end
        end

        if all_ontimes>(percent_ON*total_experiment_time) % include cell if total ontime greater than set threshold
            ON_cells=ON_cells+1;
            ON_cells_index(n)=n;
        end
end

disp (ON_cells)
disp('ON cells found')
% now use only those ON cells for rest of the analysis
ON_cells_index=ON_cells_index(~isnan(ON_cells_index));

bintime2=(30*60); %bin time in seconds, for calculating FR

% x axis for graphs of FR over whole experiment period
experiment_time_start = 0; % only change this if this experiment started late and there are less than 3 days of baseline
experiment_time_end = 3600*24*expt_end;
firstpoint=experiment_time_start+(bintime2/2);
lastpoint=experiment_time_end+(bintime2/2);
xaxis=firstpoint:bintime2:lastpoint;
xaxisdays=xaxis/(3600*24); % xaxis in days not seconds

% initialize FR matrices for all cells
all_FR_normalized_RSU=NaN(length(ON_cells_index),length(xaxis));
RSU_count=0;
FR_normalized_RSU_control=NaN(length(ON_cells_index), length(xaxis));
RSU_control_count=0;

ALL_baseline_FR = NaN(length(ON_cells_index),1); % baseline value that each cell's FR is normalized to
ALL_baseline_FR_control = NaN(length(ON_cells_index),1); % baseline value that each cell's FR is normalized to

Baseline_FR_slopes = [];
% loop through all continuous 'on' cells
RSU_number=0; % reference for how many RSU's are found and analyzed
r=0;

for m=1:length(ON_cells_index)

    n=ON_cells_index(m); % index for finding analyzed cells in original MDCELLS structure

    np_thresh=.35; % new/recent data sets
  
    if MDCELLS(n).neg_pos_time>=np_thresh 

        RSU_number=RSU_number+1;

        sorted_time=sort(MDCELLS(n).time); % just in case
         
        % filter out all spikes from 'off' times
        onTimes=sort(MDCELLS(n).onTime);
        offTimes=sort(MDCELLS(n).offTime);
        
        sorted_time(sorted_time<onTimes(1))=NaN;
        sorted_time(sorted_time>offTimes(end))=NaN;
        
        if length(onTimes)>1 
            for i=2:length(onTimes)
                current_working_time = sorted_time(sorted_time>offTimes(i-1));
                current_working_time = current_working_time(current_working_time<onTimes(i));
                start_next_off = find(sorted_time >= current_working_time(1), 1);
                end_next_off = find(sorted_time <= current_working_time(end), 1);
                sorted_time(start_next_off:end_next_off) = NaN;
                clear current_working_time
            end
        end
        % spiketimes to use for analysis
        ontime_spikes=sorted_time(~isnan(sorted_time));
        
        all_onoffTimes = sort([onTimes; offTimes]);

       % ok now moving on to plotting FR dynamics for each cell
       % normalize FR's to baseline (before MD)
           
                MYbaseline_start = 1.9*3600*24; %%%  THIS WAS 2 IN MOST OF MY ANALYSES
                MYbaseline_end = 2.4*3600*24;
              
                BL_spikes = ontime_spikes(ontime_spikes>=MYbaseline_start);
                BL_spikes = BL_spikes(BL_spikes<=MYbaseline_end);
                
                onTimes_BL = onTimes(onTimes<MYbaseline_end);
                onTimes_BL = onTimes_BL(onTimes_BL>MYbaseline_start);
                offTimes_BL = offTimes(offTimes<MYbaseline_end);
                offTimes_BL = offTimes_BL(offTimes_BL>MYbaseline_start);
 
                before_period_onoff = all_onoffTimes(all_onoffTimes<MYbaseline_start);
                after_period_onoff = all_onoffTimes(all_onoffTimes>MYbaseline_end);
                before_off = 0;
                if ~isempty(before_period_onoff)
                    if ismember(before_period_onoff(end), onTimes)
                        onTimes_BL = [MYbaseline_start; onTimes_BL];
                    elseif ismember(before_period_onoff(end), offTimes)
                        before_off = 1;
                        if isempty(onTimes_BL)   
                            disp ('This cell not on during normalization baseline!')
                            disp (n) 
                            continue
                        end
                    else
                        keyboard
                    end
                end
                if ~isempty(after_period_onoff)
                    if ismember(after_period_onoff(1), offTimes)
                        offTimes_BL = [offTimes_BL; MYbaseline_end];
                    end
                end
                if length(onTimes_BL) ~= length(offTimes_BL)
                    keyboard
                end
                if isempty(offTimes_BL)
                    continue
                end
                onTime_BL = offTimes_BL(1)-onTimes_BL(1); 
                if length(onTimes_BL)>1 
                    disp ('Length of on/off times > 1!')
                    disp (n)
                    for i=2:length(onTimes_BL)
                        onTime_BL = onTime_BL+(offTimes_BL(i)-onTimes_BL(i)); 
                    end
                end

                if onTime_BL < 3600*5 
                    disp ('This cell normalization baseline too short! SKIP')
                    disp (n)    
                    continue
                end
                baseline_FR = length(BL_spikes)/onTime_BL;
              
            % find normalized FR normalized to 'baseline' for whole
            % experiment
            edges2=experiment_time_start:bintime2:experiment_time_end;
            [bincounts2, binedges2]=histcounts(ontime_spikes, edges2);
            FR_normalized=(bincounts2/bintime2)/baseline_FR;
            
            MD_start_bintime2=round((3600*24*2.5)/bintime2); %start analyzing MD FR's compared to baseline FR; in bin number
            MD_end_bintime2=MD_start_bintime2+((3600*24*4)/bintime2); %end analyzing MD FR's for decrease 4 days after MD start
            
            % take out normalized FR during off times - put NaN's instead
            % of zeros (because already took out offtime spikes above)
            FR_normalized(1:round(onTimes(1)/bintime2))=NaN;
            FR_normalized(round(offTimes(end)/bintime2):end)=NaN;
            if length(onTimes)>1
            for ff=2:length(onTimes)
                    FR_normalized(round(offTimes(ff-1)/bintime2):round(onTimes(ff)/bintime2))=NaN;
            end
            end
            
            FR_normalized(FR_normalized>12)=NaN; % exclude normalized FR above 12
            
          
            if MDCELLS(n).deprived==1 %if cell from deprived hemisphere
                r=r+1;
                RSU_count=RSU_count+1;
                all_FR_normalized_RSU(m, 1:length(FR_normalized))= FR_normalized;
                % store baseline_FR for each cell to transform back to
                % individual FR's for ladder plots etc
                ALL_baseline_FR(m) = baseline_FR;
       
            else  %cell from non-deprived hemisphere
                r=r+1;
                RSU_control_count=RSU_control_count+1;
                FR_normalized_RSU_control(m, 1:length(FR_normalized))= FR_normalized; 
                % store baseline_FR for each cell to transform back to
                % individual FR's for ladder plots etc
                ALL_baseline_FR_control(m) = baseline_FR;             
            end
    end
end
%%


FR_means_RSU=nanmean(all_FR_normalized_RSU);
if ~isempty(FR_normalized_RSU_control(~isnan(FR_normalized_RSU_control)))
FR_means_RSU_control=nanmean(FR_normalized_RSU_control);
end

FR_SEMs_RSU=(nanstd(all_FR_normalized_RSU))./(sqrt(RSU_count)-1);
if ~isempty(FR_normalized_RSU_control(~isnan(FR_normalized_RSU_control)))
FR_SEMs_RSU_control= (nanstd(FR_normalized_RSU_control))./(sqrt(RSU_control_count)-1);
end

disp('RSU_count')
disp(RSU_count)
disp('RSU_control_count')
disp(RSU_control_count)

totalcells=RSU_count+RSU_control_count;

disp ('total cells plotted:' )
disp(totalcells)

if plot_on
figure ()
hold on
y_lims = ([0,3]);
LD_transition = 0:.5:9;
for j = 1:length(LD_transition)
    if rem(j,2)
        rectangle('Position', [LD_transition(j) 0 .5 y_lims(2)], 'FaceColor', [1 1 0 0.05])
    else
        rectangle('Position', [LD_transition(j) 0 .5 y_lims(2)], 'FaceColor', [0 0 0 0.05])
    end
end

if length(xaxisdays) > length(FR_means_RSU)
  
if ~isempty(FR_normalized_RSU_control(~isnan(FR_normalized_RSU_control)))
RSU_control=plot(xaxisdays(1:length(FR_means_RSU_control)), FR_means_RSU_control, 'k', 'DisplayName', 'control');
shadedErrorBar_JBmanuscript(xaxisdays(1:length(FR_means_RSU_control)), FR_means_RSU_control, FR_SEMs_RSU_control, 'k');
end
hold on
    
handle = plot(xaxisdays(1:length(FR_means_RSU)), FR_means_RSU, 'm', 'DisplayName', 'deprived');
shadedErrorBar_JBmanuscript(xaxisdays(1:length(FR_means_RSU)), FR_means_RSU, FR_SEMs_RSU, 'm');
else

FR_means_RSU_short = FR_means_RSU(1:length(xaxisdays));
FR_SEMs_RSU_short = FR_SEMs_RSU(1:length(xaxisdays));
    
if ~isempty(FR_normalized_RSU_control(~isnan(FR_normalized_RSU_control)))
FR_means_RSU_control_short = FR_means_RSU_control(1:length(xaxisdays));
FR_SEMs_RSU_control_short = FR_SEMs_RSU_control(1:length(xaxisdays));
RSU_control=plot(xaxisdays, FR_means_RSU_control_short, 'k', 'DisplayName', 'control');
shadedErrorBar_JBmanuscript(xaxisdays, FR_means_RSU_control_short, FR_SEMs_RSU_control_short, 'k');
end

hold on
handle = plot(xaxisdays, FR_means_RSU_short, 'm', 'DisplayName', 'deprived');
shadedErrorBar_JBmanuscript(xaxisdays, FR_means_RSU_short, FR_SEMs_RSU_short, 'm');
end

plot ([0,xaxisdays(end)], [1,1], 'k--')

ylim(y_lims);
xlim([0,xaxisdays(end)])
xlabel('Time (days)')
ylabel('Normalized Firing Rate')
%title('RSU''s: deprived vs. control')
set (gca, 'fontsize', 15)
set (gca, 'box', 'off')
end

% LADDER PLOTS V2: USING NORMALIZED FR (SINCE THIS HAS BEEN MOST THOROUGHLY
% FILTERED FOR OUTLIERS ETC, THEN TRANSFORMED BACK TO RAW FR WITH EACH
% CELL'S BASELINE VALUE
period_v2 = 12*3600; 
MD2_end_v2 = 108;
MD4_end_v2 = 168;
MD2_bins = [round((3600*MD2_end_v2)/bintime2) - round(period_v2/bintime2), round((3600*MD2_end_v2)/bintime2)];
MD4_bins = [round((3600*MD4_end_v2)/bintime2) - round(period_v2/bintime2), round((3600*MD4_end_v2)/bintime2)];

MD_ladder_deprived = NaN(length(ON_cells_index), 3);
MD_ladder_deprived(:,1) = ALL_baseline_FR;
MD_ladder_deprived(:,2) = nanmean(all_FR_normalized_RSU(:,MD2_bins(1):MD2_bins(2)),2);
MD_ladder_deprived(:,2) = MD_ladder_deprived(:,2) .* ALL_baseline_FR;
MD_ladder_deprived(:,3) = nanmean(all_FR_normalized_RSU(:,MD4_bins(1):MD4_bins(2)),2);
MD_ladder_deprived(:,3) = MD_ladder_deprived(:,3) .* ALL_baseline_FR;

MD_ladder_control = NaN(length(ON_cells_index), 3);
if ~isempty(FR_normalized_RSU_control(~isnan(FR_normalized_RSU_control)))
MD_ladder_control(:,1) = ALL_baseline_FR_control;
MD_ladder_control(:,2) = nanmean(FR_normalized_RSU_control(:,MD2_bins(1):MD2_bins(2)),2);
MD_ladder_control(:,2) = MD_ladder_control(:,2) .* ALL_baseline_FR_control;
MD_ladder_control(:,3) = nanmean(FR_normalized_RSU_control(:,MD4_bins(1):MD4_bins(2)),2);
MD_ladder_control(:,3) = MD_ladder_control(:,3) .* ALL_baseline_FR_control;
end


MD_ladder_deprived_clean = [];
MD_ladder_control_clean = [];
for mm = 1:length(MD_ladder_deprived(:,1))
    if isnan(MD_ladder_deprived(mm,1)) || isnan(MD_ladder_deprived(mm,2)) || isnan(MD_ladder_deprived(mm,3))
        continue
    else
        MD_ladder_deprived_clean(end+1,:) = MD_ladder_deprived(mm,:);
    end
    
end
for mm = 1:length(MD_ladder_control(:,1))
    if isnan(MD_ladder_control(mm,1)) || isnan(MD_ladder_control(mm,2)) || isnan(MD_ladder_control(mm,3))
        continue
    else
        MD_ladder_control_clean(end+1,:) = MD_ladder_control(mm,:);
    end
    
end

names = {'Baseline', 'MD2', 'MD4'};
x = 1:length(names);

if plot_on
figure()
for u = 1:length(MD_ladder_deprived_clean(:,1))
plot(x, MD_ladder_deprived_clean(u,:),'-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'm')
hold on
end
means = [nanmean(MD_ladder_deprived_clean(:,1)), nanmean(MD_ladder_deprived_clean(:,2)), nanmean(MD_ladder_deprived_clean(:,3))];
plot(x, means, '-md', 'MarkerSize', 20, 'MarkerFaceColor', 'm',  'LineWidth', 3)

xticks([1,2,3])
set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 15)
set(gca,'yscale','log');
ylim([10e-3, 10e1])

ylabel('Firing Rate (Spikes/Second)')
title ('Deprived Cells - CONVERTED FROM NORM')
xlim([0,4])

if ~isempty(MD_ladder_control_clean)
figure()
for u = 1:length(MD_ladder_control_clean(:,1))
plot(x, MD_ladder_control_clean(u,:), '-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3)
hold on
end
means = [nanmean(MD_ladder_control_clean(:,1)), nanmean(MD_ladder_control_clean(:,2)), nanmean(MD_ladder_control_clean(:,3))];
plot(x, means, '-kd', 'MarkerSize', 20, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 3)

xticks([1,2,3])
set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 15)
set(gca,'yscale','log');

ylabel('Firing Rate (Spikes/Second)')
title ('Control Cells - CONVERTED FROM NORM')
xlim([0,4])
end
end

CELL_DATA = struct;
CELL_DATA.deprived_ladder = MD_ladder_deprived_clean;
CELL_DATA.control_ladder = MD_ladder_control_clean;
CELL_DATA.ladder_period = period_v2;
CELL_DATA.ladder_MD2end = MD2_end_v2;
CELL_DATA.ladder_MD4end = MD4_end_v2;
CELL_DATA.bintime = bintime2;
CELL_DATA.FR_deprived = all_FR_normalized_RSU;
CELL_DATA.FR_control = FR_normalized_RSU_control;



end