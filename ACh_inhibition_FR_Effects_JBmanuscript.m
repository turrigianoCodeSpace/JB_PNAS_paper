%%% Plot FR over time for MD + Hm4di in vivo experiments, and analyze FR in control (and deprived) cells after CNO injections vs before%%%
%%% This script should make an argument for whether ACh inhibition changes
%%% FR on its own!!! 
% ACh manuscript Figure 6C

% 11/2022

%Assume:
%BL1, BL2, BL3 (*MD @ 7:30pm*), MD1, MD2, MD3, MD4...
%CNO: MD3&4
%7:30/7:30 L/D schedule

clearvars -except MDCELLS CNO_MASTER *CELL* *LFP* statetimes CAMPRAT*

Hm4di_CNO_MDCELLS = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/TlabDocs/ACh_Manuscript/PaperFigs_drafts_2022/JB_MASTER_Hm4di_CNO_CELLS_v5.mat');
MDCELLS = Hm4di_CNO_MDCELLS.MDCELLS;

%% filter to only get quality 1,2 cells
all_quals = [MDCELLS.quality];

quality_1_2=find(all_quals<3); 
MDCELLS = MDCELLS(quality_1_2);

%%

%Pick out cells that are active for some percent of total recording time to
%find continuous cells
percent_ON = .7;

expt_end=7; % how many days did the experiment last

CNO_injections = [3600*24*5, 3600*24*5.5, 3600*24*6, 3600*24*6.5];

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

% initialize arrays to hold CNO injection FR data and control data to
% compare it to
inject1_deprived = NaN(length(ON_cells_index),1);
inject2_deprived = NaN(length(ON_cells_index),1);
inject3_deprived = NaN(length(ON_cells_index),1);
inject4_deprived = NaN(length(ON_cells_index),1);
Avrg_injects_deprived = NaN(length(ON_cells_index),1);

inject1_control = NaN(length(ON_cells_index),1);
inject2_control = NaN(length(ON_cells_index),1);
inject3_control = NaN(length(ON_cells_index),1);
inject4_control = NaN(length(ON_cells_index),1);
Avrg_injects_control = NaN(length(ON_cells_index),1);

noninject1_deprived = NaN(length(ON_cells_index),1); % MD1 7:30 am
noninject2_deprived = NaN(length(ON_cells_index),1); % MD1 7:30 pm
noninject3_deprived = NaN(length(ON_cells_index),1); % MD2 7:30 am
noninject4_deprived = NaN(length(ON_cells_index),1); % MD2 7:30 pm

noninject1_control = NaN(length(ON_cells_index),1);
noninject2_control = NaN(length(ON_cells_index),1);
noninject3_control = NaN(length(ON_cells_index),1);
noninject4_control = NaN(length(ON_cells_index),1);

Beforeinject1_deprived = NaN(length(ON_cells_index),1); 
Beforeinject2_deprived = NaN(length(ON_cells_index),1); 
Beforeinject3_deprived = NaN(length(ON_cells_index),1); 
Beforeinject4_deprived = NaN(length(ON_cells_index),1); 
Avrg_Beforeinjects_deprived = NaN(length(ON_cells_index),1);

Beforeinject1_control = NaN(length(ON_cells_index),1);
Beforeinject2_control = NaN(length(ON_cells_index),1);
Beforeinject3_control = NaN(length(ON_cells_index),1);
Beforeinject4_control = NaN(length(ON_cells_index),1);
Avrg_Beforeinjects_control = NaN(length(ON_cells_index),1);

Baseline_FR_slopes = [];
% loop through all continuous 'on' cells
RSU_number=0; % reference for how many RSU's are found and analyzed
r=0;

for m=1:length(ON_cells_index)

    n=ON_cells_index(m); % index for finding analyzed cells in original MDCELLS structure
    
        np_thresh=.35; % threshold for regular spiking units
  
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
  
            %%%%%% compare firing rates before and after CNO injections %%
            %%%%%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            CNOperiod = 2*3600;
            inject1_FR = FR_normalized(round(CNO_injections(1)/bintime2):round((CNO_injections(1)+CNOperiod)/bintime2)); % light period
            inject2_FR = FR_normalized(round(CNO_injections(2)/bintime2):round((CNO_injections(2)+CNOperiod)/bintime2)); % dark period
            inject3_FR = FR_normalized(round(CNO_injections(3)/bintime2):round((CNO_injections(3)+CNOperiod)/bintime2)); % light period
            inject4_FR = FR_normalized(round(CNO_injections(4)/bintime2):round((CNO_injections(4)+CNOperiod)/bintime2)); % dark period
                     
            Beforeinject1_FR = FR_normalized(round((CNO_injections(1)-CNOperiod)/bintime2):round(CNO_injections(1)/bintime2)); % light period
            Beforeinject2_FR = FR_normalized(round((CNO_injections(2)-CNOperiod)/bintime2):round(CNO_injections(2)/bintime2)); % dark period
            Beforeinject3_FR = FR_normalized(round((CNO_injections(3)-CNOperiod)/bintime2):round(CNO_injections(3)/bintime2)); % light period
            Beforeinject4_FR = FR_normalized(round((CNO_injections(4)-CNOperiod)/bintime2):round(CNO_injections(4)/bintime2)); % dark period
            
            if MDCELLS(n).deprived==1 %if cell from deprived hemisphere
                r=r+1;
                RSU_count=RSU_count+1;
                all_FR_normalized_RSU(m, 1:length(FR_normalized))= FR_normalized;
                % store baseline_FR for each cell to transform back to
                % individual FR's for ladder plots etc
                ALL_baseline_FR(m) = baseline_FR;
                inject1_deprived(m) = nanmean(inject1_FR);
                inject2_deprived(m) = nanmean(inject2_FR);
                inject3_deprived(m) = nanmean(inject3_FR);
                inject4_deprived(m) = nanmean(inject4_FR);
                Avrg_injects_deprived(m) = nanmean([inject1_deprived(m), inject2_deprived(m), inject3_deprived(m), inject4_deprived(m)]);
                
                Beforeinject1_deprived(m) = nanmean(Beforeinject1_FR);
                Beforeinject2_deprived(m) = nanmean(Beforeinject2_FR);
                Beforeinject3_deprived(m) = nanmean(Beforeinject3_FR);
                Beforeinject4_deprived(m) = nanmean(Beforeinject4_FR);
                Avrg_Beforeinjects_deprived(m) = nanmean([Beforeinject1_deprived(m), Beforeinject2_deprived(m), Beforeinject3_deprived(m), Beforeinject4_deprived(m)]);
        
            else  %cell from non-deprived hemisphere
                r=r+1;
                RSU_control_count=RSU_control_count+1;
                FR_normalized_RSU_control(m, 1:length(FR_normalized))= FR_normalized; 
                % store baseline_FR for each cell to transform back to
                % individual FR's for ladder plots etc
                ALL_baseline_FR_control(m) = baseline_FR;
                inject1_control(m) = nanmean(inject1_FR);
                inject2_control(m) = nanmean(inject2_FR);
                inject3_control(m) = nanmean(inject3_FR);
                inject4_control(m) = nanmean(inject4_FR);
                Avrg_injects_control(m) = nanmean([inject1_control(m), inject2_control(m), inject3_control(m), inject4_control(m)]);
                
                Beforeinject1_control(m) = nanmean(Beforeinject1_FR);
                Beforeinject2_control(m) = nanmean(Beforeinject2_FR);
                Beforeinject3_control(m) = nanmean(Beforeinject3_FR);
                Beforeinject4_control(m) = nanmean(Beforeinject4_FR);
                Avrg_Beforeinjects_control(m) = nanmean([Beforeinject1_control(m), Beforeinject2_control(m), Beforeinject3_control(m), Beforeinject4_control(m)]);
              
            end
    end
end
%%

% COMPARE TIME RIGHT BEFORE TO RIGHT AFTER CNO INJECTIONS

BeforeAfter_inject1_control = [];
BeforeAfter_inject2_control = [];
BeforeAfter_inject3_control = [];
BeforeAfter_inject4_control = [];
BeforeAfter_avrgInjects_control = [];

BeforeAfter_inject1_deprived = [];
BeforeAfter_inject2_deprived = [];
BeforeAfter_inject3_deprived = [];
BeforeAfter_inject4_deprived = [];
BeforeAfter_avrgInjects_deprived = [];

for i = 1:length(ON_cells_index)
    if ~isnan(Beforeinject1_control(i))
        if ~isnan(inject1_control(i))
            BeforeAfter_inject1_control(end+1,:) = [Beforeinject1_control(i), inject1_control(i)];
        end
    end
    if ~isnan(Beforeinject2_control(i))
        if ~isnan(inject2_control(i))
            BeforeAfter_inject2_control(end+1,:) = [Beforeinject2_control(i), inject2_control(i)];
        end
    end
    if ~isnan(Beforeinject3_control(i))
        if ~isnan(inject3_control(i))
            BeforeAfter_inject3_control(end+1,:) = [Beforeinject3_control(i), inject3_control(i)];
        end
    end
    if ~isnan(Beforeinject4_control(i))
        if ~isnan(inject4_control(i))
            BeforeAfter_inject4_control(end+1,:) = [Beforeinject4_control(i), inject4_control(i)];
        end    
    end
    
    if ~isnan(Beforeinject1_deprived(i))
        if ~isnan(inject1_deprived(i))
            BeforeAfter_inject1_deprived(end+1,:) = [Beforeinject1_deprived(i), inject1_deprived(i)];
        end
    end
    if ~isnan(Beforeinject2_deprived(i))
        if ~isnan(inject2_deprived(i))
            BeforeAfter_inject2_deprived(end+1,:) = [Beforeinject2_deprived(i), inject2_deprived(i)];
        end
    end
    if ~isnan(Beforeinject3_deprived(i))
        if ~isnan(inject3_deprived(i))
            BeforeAfter_inject3_deprived(end+1,:) = [Beforeinject3_deprived(i), inject3_deprived(i)];
        end
    end
    if ~isnan(Beforeinject4_deprived(i))
        if ~isnan(inject4_deprived(i))
            BeforeAfter_inject4_deprived(end+1,:) = [Beforeinject4_deprived(i), inject4_deprived(i)];
        end    
    end
    if ~isnan(Avrg_Beforeinjects_deprived(i))
        if ~isnan(Avrg_injects_deprived(i))
            BeforeAfter_avrgInjects_deprived(end+1,:) = [Avrg_Beforeinjects_deprived(i), Avrg_injects_deprived(i)];
        end    
    end
    if ~isnan(Avrg_Beforeinjects_control(i))
        if ~isnan(Avrg_injects_control(i))
            BeforeAfter_avrgInjects_control(end+1,:) = [Avrg_Beforeinjects_control(i), Avrg_injects_control(i)];
        end    
    end
    
end

% CONTROL HEM CELLS
figure()
hold on

for u = 1:length(BeforeAfter_inject1_control(:,1))
plot([1,2], BeforeAfter_inject1_control(u,:), '-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3)
hold on
end
for u = 1:length(BeforeAfter_inject2_control(:,1))
plot([1,2], BeforeAfter_inject2_control(u,:), '-bo', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3)
hold on
end
for u = 1:length(BeforeAfter_inject3_control(:,1))
plot([1,2], BeforeAfter_inject3_control(u,:), '-go', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3)
hold on
end
for u = 1:length(BeforeAfter_inject4_control(:,1))
plot([1,2], BeforeAfter_inject4_control(u,:), '-co', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3)
hold on
end
means1 = [nanmean(BeforeAfter_inject1_control(:,1)), nanmean(BeforeAfter_inject1_control(:,2))];
plot([1,2], means1, '-kd', 'MarkerSize', 20, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 3)
means2 = [nanmean(BeforeAfter_inject2_control(:,1)), nanmean(BeforeAfter_inject2_control(:,2))];
plot([1,2], means2, '-bd', 'MarkerSize', 20, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'LineWidth', 3)
means3 = [nanmean(BeforeAfter_inject3_control(:,1)), nanmean(BeforeAfter_inject3_control(:,2))];
plot([1,2], means3, '-gd', 'MarkerSize', 20, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'LineWidth', 3)
means4 = [nanmean(BeforeAfter_inject4_control(:,1)), nanmean(BeforeAfter_inject4_control(:,2))];
plot([1,2], means4, '-cd', 'MarkerSize', 20, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'c', 'LineWidth', 3)

ALL_means_BeforeAfter = [BeforeAfter_inject1_control;  BeforeAfter_inject2_control; BeforeAfter_inject3_control; BeforeAfter_inject4_control];
ALL_means = [nanmean(ALL_means_BeforeAfter(:,1)), nanmean(ALL_means_BeforeAfter(:,2))];
plot([1,2], ALL_means, '-rd', 'MarkerSize', 20, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineWidth', 3)
[p_control,h] = signrank(ALL_means_BeforeAfter(:,1), ALL_means_BeforeAfter(:,2))

xticks([1,2])
set (gca, 'xticklabel', {'Before CNO inj', 'After CNO inj'}, 'box', 'off', 'fontsize', 15)
set(gca,'yscale','log');

ylabel('Firing Rate (Spikes/Second)')
title ('Before/After: Control Cells')
xlim([0,3])

% DEPRIVED HEM CELLS
figure()
hold on

for u = 1:length(BeforeAfter_inject1_deprived(:,1))
plot([1,2], BeforeAfter_inject1_deprived(u,:), '-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3)
hold on
end
for u = 1:length(BeforeAfter_inject2_deprived(:,1))
plot([1,2], BeforeAfter_inject2_deprived(u,:), '-bo', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3)
hold on
end
for u = 1:length(BeforeAfter_inject3_deprived(:,1))
plot([1,2], BeforeAfter_inject3_deprived(u,:), '-go', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3)
hold on
end
for u = 1:length(BeforeAfter_inject4_deprived(:,1))
plot([1,2], BeforeAfter_inject4_deprived(u,:), '-co', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3)
hold on
end
means1 = [nanmean(BeforeAfter_inject1_deprived(:,1)), nanmean(BeforeAfter_inject1_deprived(:,2))];
plot([1,2], means1, '-kd', 'MarkerSize', 20, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 3)
means2 = [nanmean(BeforeAfter_inject2_deprived(:,1)), nanmean(BeforeAfter_inject2_deprived(:,2))];
plot([1,2], means2, '-bd', 'MarkerSize', 20, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'LineWidth', 3)
means3 = [nanmean(BeforeAfter_inject3_deprived(:,1)), nanmean(BeforeAfter_inject3_deprived(:,2))];
plot([1,2], means3, '-gd', 'MarkerSize', 20, 'MarkerFaceColor', 'g', 'MarkerEdgeColor', 'g', 'LineWidth', 3)
means4 = [nanmean(BeforeAfter_inject4_deprived(:,1)), nanmean(BeforeAfter_inject4_deprived(:,2))];
plot([1,2], means4, '-cd', 'MarkerSize', 20, 'MarkerFaceColor', 'c', 'MarkerEdgeColor', 'c', 'LineWidth', 3)

ALL_means_BeforeAfter = [BeforeAfter_inject1_deprived;  BeforeAfter_inject2_deprived; BeforeAfter_inject3_deprived; BeforeAfter_inject4_deprived];
ALL_means = [nanmean(ALL_means_BeforeAfter(:,1)), nanmean(ALL_means_BeforeAfter(:,2))];
plot([1,2], ALL_means, '-rd', 'MarkerSize', 20, 'MarkerFaceColor', 'r', 'MarkerEdgeColor', 'r', 'LineWidth', 3)
[p_deprived,h] = signrank(ALL_means_BeforeAfter(:,1), ALL_means_BeforeAfter(:,2))

xticks([1,2])
set (gca, 'xticklabel', {'Before CNO inj', 'After CNO inj'}, 'box', 'off', 'fontsize', 15)
set(gca,'yscale','log');

ylabel('Firing Rate (Spikes/Second)')
title ('Before/After: Deprived Cells')
xlim([0,3])

%%%% now before / after with average of all injections per cell
% CONTROL HEM CELLS
figure()
hold on

for u = 1:length(BeforeAfter_avrgInjects_control(:,1))
plot([1,2], BeforeAfter_avrgInjects_control(u,:), '-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3)
hold on
end

ALL_means = [nanmean(BeforeAfter_avrgInjects_control(:,1)), nanmean(BeforeAfter_avrgInjects_control(:,2))];
plot([1,2], ALL_means, '-kd', 'MarkerSize', 20, 'MarkerFaceColor', 'k', 'MarkerEdgeColor', 'k', 'LineWidth', 3)
[p_controlALL,h] = signrank(BeforeAfter_avrgInjects_control(:,1), BeforeAfter_avrgInjects_control(:,2))

xticks([1,2])
set (gca, 'xticklabel', {'Before CNO inj', 'After CNO inj'}, 'box', 'off', 'fontsize', 15)
set(gca,'yscale','log');

ylabel('Firing Rate (Spikes/Second)')
title ('Before/After: Control Cells (AVRG per cell)')
xlim([0,3])

% DEPRIVED HEM CELLS
figure()
hold on

for u = 1:length(BeforeAfter_avrgInjects_deprived(:,1))
plot([1,2], BeforeAfter_avrgInjects_deprived(u,:), '-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'MarkerEdgeColor', 'b', 'LineWidth', 3)
hold on
end

ALL_means = [nanmean(BeforeAfter_avrgInjects_deprived(:,1)), nanmean(BeforeAfter_avrgInjects_deprived(:,2))];
plot([1,2], ALL_means, '-kd', 'MarkerSize', 20, 'MarkerFaceColor', 'b', 'MarkerEdgeColor', 'b', 'LineWidth', 3)

[p_deprivedALL,h] = signrank(BeforeAfter_avrgInjects_deprived(:,1), BeforeAfter_avrgInjects_deprived(:,2))

xticks([1,2])
set (gca, 'xticklabel', {'Before CNO inj', 'After CNO inj'}, 'box', 'off', 'fontsize', 15)
set(gca,'yscale','log');

ylabel('Firing Rate (Spikes/Second)')
title ('Before/After: Deprived Cells (AVRG per cell)')
xlim([0,3])
ylim([10e-3,10e0])
