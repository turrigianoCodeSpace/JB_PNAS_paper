%%% Plot EXAMPLE FR TRACES FOR HM4DI+CNO AND CNO ONLY GROUPS %%%
 
%Assume:
%BL1, BL2, BL3 (*MD @ 7:30pm*), MD1, MD2, MD3, MD4...
%if CNO: MD3&4
%7:30/7:30 L/D schedule

clearvars -except MDCELLS CNO_MASTER *CELL* *LFP* statetimes 


%%

Hm4di_cell = 83;
CNOonly_cell = 37;
 
bintime2=(30*60); %bin time in seconds, for calculating FR

% x axis for graphs of FR over whole experiment period
experiment_time_start = 0; % only change this if this experiment started late and there are less than 3 days of baseline
experiment_time_end = 3600*24*7;
firstpoint=experiment_time_start+(bintime2/2);
lastpoint=experiment_time_end+(bintime2/2);
xaxis=firstpoint:bintime2:lastpoint;
xaxisdays=xaxis/(3600*24); % xaxis in days not seconds


for m=1:2
    clear MDCELLS
    if m == 1
        Hm4di_MDCELLS = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/TlabDocs/ACh_Manuscript/PaperFigs_drafts_2022/JB_MASTER_Hm4di_CNO_CELLS_v5.mat');
        MDCELLS = Hm4di_MDCELLS.MDCELLS;
        n = Hm4di_cell;
        this_color = [0,0,1];
    elseif m == 2
        CNOonly_MDCELLS = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/TlabDocs/ACh_Manuscript/PaperFigs_drafts_2022/JB_MASTER_CNOonly_CELLS_v2.mat');
        MDCELLS = CNOonly_MDCELLS.MDCELLS;
        n = CNOonly_cell;
        this_color = [1,0,1];
    end
    
    % to make sure cell numbering is consistent with other scripts
    all_quals = [MDCELLS.quality];
    quality_1_2=find(all_quals<3); 
    MDCELLS = MDCELLS(quality_1_2);

edges2=experiment_time_start:bintime2:experiment_time_end;
[bincounts2, binedges2]=histcounts(MDCELLS(n).time, edges2);
this_FR=bincounts2/bintime2;

% for plotting baseline reference line
baseline_start_bintime2 = round((3600*24*2)/bintime2);
baseline_end_bintime2 = round((3600*24*2.4)/bintime2);

baseline_FR = mean(this_FR(baseline_start_bintime2:baseline_end_bintime2));

figure(1)
if length(xaxisdays) > length(this_FR)
    plot(xaxisdays(1:length(this_FR)), this_FR, 'Color', this_color, 'LineWidth', 2)
else 
    this_FR_short = this_FR(1:length(xaxisdays));
    plot(xaxisdays, this_FR_short, 'Color', this_color, 'LineWidth', 2)
end
hold on 

plot([experiment_time_start, experiment_time_end], [baseline_FR, baseline_FR], 'k--')

xlabel('Time (days)')
ylabel('Firing Rate (Hz)')
set (gca, 'fontsize', 15, 'box', 'off', 'xlim', [2,7])
set(gca,'yscale','log');

y_lims = ([10e-2,10e1]);
LD_transition = 0:.5:9;
for j = 1:length(LD_transition)
    if rem(j,2)
        rectangle('Position', [LD_transition(j) y_lims(1) .5 y_lims(2)], 'FaceColor', [1 1 0 0.05])
    else
        rectangle('Position', [LD_transition(j) y_lims(1) .5 y_lims(2)], 'FaceColor', [0 0 0 0.05])
    end
end
ylim(y_lims)
            
end



