% plot aspects of depolarizations and hyperpolarizations from CCh response
% experiments
% ACh manuscript - CompareCChResponses_JBmanuscript.m calls this script for
% Figure 4F

function LADDERS = CChResponse_Plots_JBmanuscript(condition, plot_on)

if nargin == 1
    plot_on = 1;
end

if strcmp(condition, 'ctrl') || strcmp(condition, 'CTRL')
    DATA = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/ALL_ctrl_CCh_Analysis.mat');
elseif strcmp(condition, 'm1_3') || strcmp(condition, 'M1_3')
    DATA = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/m1_3sh_CCh_Analysis.mat');
else
    fprintf('ERROR: UNKOWN CONDITION; TRY "CTRL" or "m1_3" \n')
end

DATA = DATA.FI_DATA;

% choose sweeps to analyze/compare
baseline = [1:25]; 
CCh = [125:150];
wash = [225:250];

% initiate variables of interest
all_lengths = [];
for i = 1:length(DATA.filename)
    all_lengths = [all_lengths, length(DATA.Ra{1,i})];
    all_lengths = [all_lengths, length(DATA.Rin_IC{1,i})];
    all_lengths = [all_lengths, length(DATA.MFR{1,i})];
end

ALL_PP = struct;
ALL_PP.data = nan(length(DATA.filename), 5);
ALL_PP.labels = {'Ra_Vc', 'Rin_Vc', 'Rin_IC', 'Cp_Vc', 'Cp_IC'};

DATA_IFR_ave = nan(length(DATA.filename), max(all_lengths));
LADDER_IFR_ave = nan(length(DATA.filename), 3);

DATA_MFR = nan(length(DATA.filename), max(all_lengths));
LADDER_MFR = nan(length(DATA.filename), 3);

DATA_Rin = nan(length(DATA.filename), max(all_lengths));
DATA_Rin_VC = nan(length(DATA.filename), max(all_lengths));
LADDER_Rin = nan(length(DATA.filename), 3);

DATA_Cp = nan(length(DATA.filename), max(all_lengths));
DATA_Cp_VC = nan(length(DATA.filename), max(all_lengths));
LADDER_Cp = nan(length(DATA.filename), 3);

DATA_Ra = nan(length(DATA.filename), max(all_lengths));

DATA_latency = nan(length(DATA.filename), max(all_lengths));

DATA_width = nan(length(DATA.filename), max(all_lengths));

DATA_Vrm = nan(length(DATA.filename), max(all_lengths));
LADDER_Vrm = nan(length(DATA.filename), 3);

DATA_adpind = nan(length(DATA.filename), max(all_lengths));

DATA_AHP = nan(length(DATA.filename), max(all_lengths));
DATA_ADP = nan(length(DATA.filename), max(all_lengths));
LADDER_ADP = nan(length(DATA.filename), 3);
LADDER_AHP = nan(length(DATA.filename), 3);

NAMES = cell(1, length(DATA.filename));
% only for first two test days - flow rate 1 vs flow rate 2 (ended up
% choosing 2)
FlowRate_1 = {'date063022_cell_03', 'date063022_cell_04', 'date063022_cell_08', 'date062922_cell_02', 'date062922_cell_03', 'date062922_cell_01'}; % 'date063022_cell_01' was also flow rate 1, but lost access completely (~26MO by the end)
FlowRate_2 = {'date063022_cell_02', 'date063022_cell_05', 'date063022_cell_06', 'date063022_cell_07', 'date062922_cell_04', 'date062922_cell_05', 'date062922_cell_06', 'date062922_cell_07'};

% loop through to gather data from each cell
count = 0;
for h = 1:length(DATA.filename)
    
    
    if contains(DATA.filename{1,h}, {'date063022'}) % these cells had slightly different protocol - only 200pA pos pulse instead of 250pA
        continue
    end
    
    % exclude cells with bad PP
    this_VC_Ra = nonzeros(DATA.Ra{1,h});
    this_VC_Rin = nonzeros(DATA.Rin{1,h});
    this_VC_Cp = nonzeros(DATA.Cp{1,h});
    
    this_IC_Rin = nonzeros(DATA.Rin_IC{1,h});
    this_IC_Cp = nonzeros(DATA.Cp_IC{1,h});
    
    ALL_PP.data(h, :) = [mean(this_VC_Ra), mean(this_VC_Rin), mean(this_IC_Rin), mean(this_VC_Cp), mean(this_IC_Cp)];
    

    if ~isempty(this_VC_Ra)
        if mean(this_VC_Ra) > 20000000
            fprintf('\n Cell %s too high Ra \n', DATA.filename{1,h})
            continue
        elseif mean(this_VC_Rin) < 90000000 || mean(this_VC_Rin) > 350000000
            fprintf('\n Cell %s too low/high Rin \n', DATA.filename{1,h})
            disp (this_VC_Rin .* 10^-6)
            continue
        end
    end
    
    % check to make sure cell was recorded long enough to measure washout
    if length(DATA.Rin_IC{1,h})<wash(2)
        if length(DATA.Rin_IC{1,h})>wash(1) + 9
            this_wash = [wash(1), length(DATA.Rin_IC{1,h})];
        else
            this_wash = 0;
        end
    else
        this_wash = wash;
    end

    count = count + 1;
    
    NAMES{count} = DATA.filename{1,h};
    
    IC_sealTest_traces = find(DATA.Rin_IC{1,h} ~= 0);
    VC_sealTest_traces = find(DATA.Ra{1,h} ~= 0);
    
    Rin_IC_zeros = find(DATA.Rin_IC{1,h} == 0);
    thisRin_IC = DATA.Rin_IC{1,h};
    thisRin_IC(Rin_IC_zeros) = NaN;
    
    Cp_IC_zeros = find(DATA.Cp_IC{1,h} == 0);
    thisCp_IC = DATA.Cp_IC{1,h};
    thisCp_IC(Cp_IC_zeros) = NaN;
    
    DATA_IFR_ave(count, 1:length(DATA.IFR_ave{1,h})) = DATA.IFR_ave{1,h}';
    DATA_IFR_ave(count, IC_sealTest_traces) = NaN;
    DATA_IFR_ave(count, VC_sealTest_traces) = NaN;
    if ismember(DATA.filename{1,h}, {'date071322_cell_07', 'date071522_cell_05'})
        DATA_IFR_ave(count, 2:length(DATA.IFR_ave{1,h})-1) = DATA_IFR_ave(count, 3:length(DATA.IFR_ave{1,h}));
        DATA_IFR_ave(count, length(DATA.IFR_ave{1,h})) = NaN;
    end
    
    DATA_MFR(count, 1:length(DATA.MFR{1,h})) = DATA.MFR{1,h}';
    DATA_MFR(count, IC_sealTest_traces) = NaN;
    DATA_MFR(count, VC_sealTest_traces) = NaN;
    if ismember(DATA.filename{1,h}, {'date071322_cell_07', 'date071522_cell_05'})
        DATA_MFR(count, 2:length(DATA.MFR{1,h})-1) = DATA_MFR(count, 3:length(DATA.MFR{1,h}));
        DATA_MFR(count, length(DATA.MFR{1,h})) = NaN;
    end
    
    DATA_latency(count, 1:length(DATA.lat{1,h})) = DATA.lat{1,h}';
    DATA_latency(count, IC_sealTest_traces) = NaN;
    DATA_latency(count, VC_sealTest_traces) = NaN;
    if ismember(DATA.filename{1,h}, {'date071322_cell_07', 'date071522_cell_05'})
        DATA_latency(count, 2:length(DATA.lat{1,h})-1) = DATA_latency(count, 3:length(DATA.lat{1,h}));
        DATA_latency(count, length(DATA.lat{1,h})) = NaN;
    end
    
    DATA_width(count, 1:length(DATA.width_median{1,h})) = DATA.width_median{1,h}';
    DATA_width(count, IC_sealTest_traces) = NaN;
    DATA_width(count, VC_sealTest_traces) = NaN;
    if ismember(DATA.filename{1,h}, {'date071322_cell_07', 'date071522_cell_05'})
        DATA_width(count, 2:length(DATA.width_median{1,h})-1) = DATA_width(count, 3:length(DATA.width_median{1,h}));
        DATA_width(count, length(DATA.width_median{1,h})) = NaN;
    end
    
    DATA_Vrm(count, 1:length(DATA.Vrm{1,h})) = DATA.Vrm{1,h}';
    DATA_adpind(count, 1:length(DATA.adp_index{1,h})) = DATA.adp_index{1,h}';
    DATA_adpind(count, IC_sealTest_traces) = NaN;
    DATA_adpind(count, VC_sealTest_traces) = NaN;
    if ismember(DATA.filename{1,h}, {'date071322_cell_07', 'date071522_cell_05'})
        DATA_Vrm(count, 2:length(DATA.Vrm{1,h})-1) = DATA_Vrm(count, 3:length(DATA.Vrm{1,h}));
        DATA_Vrm(count, length(DATA.Vrm{1,h})) = NaN;
    end
    
    DATA_ADP(count, 1:length(DATA.ADP{1,h})) = DATA.ADP{1,h}';
    DATA_AHP(count, 1:length(DATA.AHP{1,h})) = DATA.AHP{1,h}';
    DATA_Rin(count, 1:length(DATA.Rin_IC{1,h})) = thisRin_IC';
    DATA_Rin_VC(count, 1:length(DATA.Rin{1,h})) = DATA.Rin{1,h}';
    DATA_Cp(count, 1:length(DATA.Cp_IC{1,h})) = thisCp_IC';
    DATA_Cp_VC(count, 1:length(DATA.Cp{1,h})) = DATA.Cp{1,h}';
    DATA_Ra(count, 1:length(DATA.Ra{1,h})) = DATA.Ra{1,h}';
    
    if ismember(DATA.filename{1,h}, {'date071322_cell_07', 'date071522_cell_05'})
        DATA_ADP(count, 2:length(DATA.ADP{1,h})-1) = DATA_ADP(count, 3:length(DATA.ADP{1,h}));
        DATA_ADP(count, length(DATA.ADP{1,h})) = NaN;
        
        DATA_AHP(count, 2:length(DATA.AHP{1,h})-1) = DATA_AHP(count, 3:length(DATA.AHP{1,h}));
        DATA_AHP(count, length(DATA.AHP{1,h})) = NaN;
        
        DATA_Rin(count, 2:length(DATA.Rin_IC{1,h})-1) = DATA_Rin(count, 3:length(DATA.Rin_IC{1,h}));
        DATA_Rin(count, length(DATA.Rin_IC{1,h})) = NaN;
        
        DATA_Rin_VC(count, 2:length(DATA.Rin{1,h})-1) = DATA_Rin_VC(count, 3:length(DATA.Rin{1,h}));
        DATA_Rin_VC(count, length(DATA.Rin{1,h})) = NaN;
        
        DATA_Cp(count, 2:length(DATA.Cp_IC{1,h})-1) = DATA_Cp(count, 3:length(DATA.Cp_IC{1,h}));
        DATA_Cp(count, length(DATA.Cp_IC{1,h})) = NaN;
        
        DATA_Cp_VC(count, 2:length(DATA.Cp{1,h})-1) = DATA_Cp_VC(count, 3:length(DATA.Cp{1,h}));
        DATA_Cp_VC(count, length(DATA.Cp{1,h})) = NaN;
    end
    
    % gather data for ladder plots
    if this_wash ~= 0
    LADDER_Rin(h, :) = [nanmean(thisRin_IC(baseline)), nanmean(thisRin_IC(CCh)), nanmean(thisRin_IC(this_wash))];
    LADDER_Cp(h, :) = [nanmean(thisCp_IC(baseline)), nanmean(thisCp_IC(CCh)), nanmean(thisCp_IC(this_wash))];
    LADDER_IFR_ave(h, :) = [nanmean(DATA_IFR_ave(count, baseline)), nanmean(DATA_IFR_ave(count, CCh)), nanmean(DATA_IFR_ave(count, this_wash))];
    LADDER_MFR(h, :) = [nanmean(DATA_MFR(count, baseline)), nanmean(DATA_MFR(count, CCh)), nanmean(DATA_MFR(count, this_wash))];
    LADDER_ADP(h, :) = [nanmean(DATA_ADP(count, baseline)), nanmean(DATA_ADP(count, CCh)), nanmean(DATA_ADP(count, this_wash))];
    LADDER_AHP(h, :) = [nanmean(DATA_AHP(count, baseline)), nanmean(DATA_AHP(count, CCh)), nanmean(DATA_AHP(count, this_wash))];
    LADDER_Vrm(h, :) = [nanmean(DATA_Vrm(count, baseline)), nanmean(DATA_Vrm(count, CCh)), nanmean(DATA_Vrm(count, this_wash))];
    else
    LADDER_Rin(h, :) = [nanmean(thisRin_IC(baseline)), nanmean(thisRin_IC(CCh)), NaN];
    LADDER_Cp(h, :) = [nanmean(thisCp_IC(baseline)), nanmean(thisCp_IC(CCh)), NaN];
    LADDER_IFR_ave(h, :) = [nanmean(DATA_IFR_ave(count, baseline)), nanmean(DATA_IFR_ave(count, CCh)), NaN];
    LADDER_MFR(h, :) = [nanmean(DATA_MFR(count, baseline)), nanmean(DATA_MFR(count, CCh)), NaN];
    LADDER_ADP(h, :) = [nanmean(DATA_ADP(count, baseline)), nanmean(DATA_ADP(count, CCh)), NaN];
    LADDER_AHP(h, :) = [nanmean(DATA_AHP(count, baseline)), nanmean(DATA_AHP(count, CCh)), NaN];
    LADDER_Vrm(h, :) = [nanmean(DATA_Vrm(count, baseline)), nanmean(DATA_Vrm(count, CCh)), NaN];
    end

end

Colors = colormap(jet(count));



mean_IFR = nanmean(DATA_IFR_ave, 1);
sem_IFR = nanstd(DATA_IFR_ave)./sqrt(count-1);

mean_MFR = nanmean(DATA_MFR, 1);
sem_MFR = nanstd(DATA_MFR)./sqrt(count-1);

mean_latency = nanmean(DATA_latency, 1);
sem_latency = nanstd(DATA_latency)./sqrt(count-1);

mean_width = nanmean(DATA_width, 1);
sem_width = nanstd(DATA_width)./sqrt(count-1);

mean_Vrm = nanmean(DATA_Vrm, 1);
sem_Vrm = nanstd(DATA_Vrm)./sqrt(count-1);

mean_adpind = nanmean(DATA_adpind, 1);
sem_adpind = nanstd(DATA_adpind)./sqrt(count-1);

mean_AHP = nanmean(DATA_AHP, 1);
sem_AHP = nanstd(DATA_AHP)./sqrt(count-1);

mean_ADP = nanmean(DATA_ADP, 1);
sem_ADP = nanstd(DATA_ADP)./sqrt(count-1);

mean_Rin = nanmean(DATA_Rin, 1);
sem_Rin = nanstd(DATA_Rin)./sqrt(count-1);

mean_Cp = nanmean(DATA_Cp, 1);
sem_Cp = nanstd(DATA_Cp)./sqrt(count-1);



% plot!

if plot_on

% Plot Rin
figure()
hold on
for f = 1:count
    plot(DATA_Rin(f,:) .* 10^-6, 'x', 'Color', Colors(f,:))
    plot(DATA_Rin_VC(f,:) .* 10^-6 , 'o', 'Color', Colors(f,:))
end
plot ([25, 25], [50, 200], 'k', 'LineWidth', 3)
plot ([150, 150], [50, 200], 'k', 'LineWidth', 3)

x_meanRin = find(~isnan(mean_Rin));
y_meanRin = ~isnan(mean_Rin);
y_meanRin = mean_Rin(y_meanRin);

x_semRin = find(~isnan(sem_Rin));
y_semRin = ~isnan(sem_Rin);
y_semRin = sem_Rin(y_semRin);
plot (x_meanRin, y_meanRin .* 10^-6, 'k', 'LineWidth', 3)
shadedErrorBar (x_meanRin, y_meanRin .* 10^-6, y_semRin .* 10^-6, 'k')

set (gca, 'box', 'off', 'fontsize', 18)
xlabel('Sweep Number')
ylabel('Rin (MOhms)')
%ylim([60,180])
ylim([60, 180])


% Plot Cp
figure()
hold on
for f = 1:count
    plot(DATA_Cp(f,:) .* 10^12, 'x', 'Color', Colors(f,:))
    plot(DATA_Cp_VC(f,:) .* 10^12 , 'o', 'Color', Colors(f,:))
end
plot ([25, 25], [50, 200], 'k', 'LineWidth', 3)
plot ([150, 150], [50, 200], 'k', 'LineWidth', 3)
plot (mean_Cp .* 10^12, 'k--', 'LineWidth', 3)

set (gca, 'box', 'off', 'fontsize', 18)
xlabel('Sweep Number')
ylabel('Cp (pF)')
ylim([60,180])

% Plot IFR_ave
figure()
hold on
for f = 1:count
    plot(DATA_IFR_ave(f,:), 'x', 'Color', Colors(f,:))
end
plot ([25, 25], [0, 30], 'k', 'LineWidth', 3)
plot ([150, 150], [0, 30], 'k', 'LineWidth', 3)

x_meanIFR = find(~isnan(mean_IFR));
y_meanIFR = ~isnan(mean_IFR);
y_meanIFR = mean_IFR(y_meanIFR);

x_semIFR = find(~isnan(sem_IFR));
y_semIFR = ~isnan(sem_IFR);
y_semIFR = sem_IFR(y_semIFR);
plot (x_meanIFR, y_meanIFR, 'k', 'LineWidth', 3)
shadedErrorBar (x_meanIFR, y_meanIFR, y_semIFR, 'k')

set (gca, 'box', 'off', 'fontsize', 18)
xlabel('Sweep Number')
ylabel('Mean IFR')
legend(NAMES{1:count})


% Plot MFR
figure()
hold on
for f = 1:count
    plot(DATA_MFR(f,:), 'x', 'Color', Colors(f,:))
end
plot ([25, 25], [0, 10], 'k', 'LineWidth', 3)
plot ([150, 150], [0, 10], 'k', 'LineWidth', 3)

x_meanMFR = find(~isnan(mean_MFR));
y_meanMFR = ~isnan(mean_MFR);
y_meanMFR = mean_MFR(y_meanMFR);

x_semMFR = find(~isnan(sem_MFR));
y_semMFR = ~isnan(sem_MFR);
y_semMFR = sem_MFR(y_semMFR);
plot (x_meanMFR, y_meanMFR, 'k', 'LineWidth', 3)
shadedErrorBar (x_meanMFR, y_meanMFR, y_semMFR, 'k')

set (gca, 'box', 'off', 'fontsize', 18)
xlabel('Sweep Number')
ylabel('Mean Firing Rate')
legend(NAMES{1:count})

% Plot Vrm
figure()
hold on
for f = 1:count
    plot(DATA_Vrm(f,:), 'x', 'Color', Colors(f,:))
end
plot ([25, 25], [-85, -60], 'k', 'LineWidth', 3)
plot ([150, 150], [-85, -60], 'k', 'LineWidth', 3)
x_meanVrm = find(~isnan(mean_Vrm));
y_meanVrm = ~isnan(mean_Vrm);
y_meanVrm = mean_Vrm(y_meanVrm);

x_semVrm = find(~isnan(sem_Vrm));
y_semVrm = ~isnan(sem_Vrm);
y_semVrm = sem_Vrm(y_semVrm);
plot (x_meanVrm, y_meanVrm, 'k', 'LineWidth', 3)
shadedErrorBar(x_meanVrm, y_meanVrm, y_semVrm, 'k')

set (gca, 'box', 'off', 'fontsize', 18)
xlabel('Sweep Number')
ylabel('Resting Membrane Potential (mV)')
legend(NAMES{1:count})


% Plot AHP
figure()
hold on
for f = 1:count
    plot(DATA_AHP(f,:), 'x', 'Color', Colors(f,:))
end
plot ([25, 25], [-5, 5], 'k', 'LineWidth', 3)
plot ([150, 150], [-5, 5], 'k', 'LineWidth', 3)
x_meanAHP = find(~isnan(mean_AHP));
y_meanAHP = ~isnan(mean_AHP);
y_meanAHP = mean_AHP(y_meanAHP);

x_semAHP = find(~isnan(sem_AHP));
y_semAHP = ~isnan(sem_AHP);
y_semAHP = sem_AHP(y_semAHP);
plot (x_meanAHP, y_meanAHP, 'k', 'LineWidth', 3)
shadedErrorBar(x_meanAHP, y_meanAHP, y_semAHP, 'k')

set (gca, 'box', 'off', 'fontsize', 18)
xlabel('Sweep Number')
ylabel('AHP (mV)')
legend(NAMES{1:count})


% Plot ADP
figure()
hold on
for f = 1:count
    plot(DATA_ADP(f,:), 'x', 'Color', Colors(f,:), 'DisplayName', NAMES{f})
end
plot ([25, 25], [-5, 5], 'k', 'LineWidth', 3)
plot ([150, 150], [-5, 5], 'k', 'LineWidth', 3)
x_meanADP = find(~isnan(mean_ADP));
y_meanADP = ~isnan(mean_ADP);
y_meanADP = mean_ADP(y_meanADP);

x_semADP = find(~isnan(sem_ADP));
y_semADP = ~isnan(sem_ADP);
y_semADP = sem_ADP(y_semADP);
plot (x_meanADP, y_meanADP, 'k', 'LineWidth', 3)
shadedErrorBar(x_meanADP, y_meanADP, y_semADP, 'k')


set (gca, 'box', 'off', 'fontsize', 18)
xlabel('Sweep Number')
ylabel('ADP (mV)')
legend(NAMES{1:count})


% plot ladder plots
names = {'Baseline', 'CCh', 'Wash'};
x = 1:length(names);

% IFR_ave LADDER

figure()
for u = 1:length(LADDER_IFR_ave(:,1))
plot(x, LADDER_IFR_ave(u,:),'-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'k')
hold on
end
means = [nanmean(LADDER_IFR_ave(:,1)), nanmean(LADDER_IFR_ave(:,2)), nanmean(LADDER_IFR_ave(:,3))];
plot(x, means, '-md', 'MarkerSize', 20, 'MarkerFaceColor', 'm',  'LineWidth', 3)

xticks([1,2,3])
set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 15)

ylabel('Mean Inst. Firing Rate')
% paired non-parametric test for baseline vs CCh - should multiply by 3
% maybe for full ladder??
p_IFR_ave = signrank(LADDER_IFR_ave(:,1), LADDER_IFR_ave(:,2))

% MFR LADDER

figure()
for u = 1:length(LADDER_MFR(:,1))
plot(x, LADDER_MFR(u,:),'-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'k')
hold on
end
means = [nanmean(LADDER_MFR(:,1)), nanmean(LADDER_MFR(:,2)), nanmean(LADDER_MFR(:,3))];
plot(x, means, '-md', 'MarkerSize', 20, 'MarkerFaceColor', 'm',  'LineWidth', 3)

xticks([1,2,3])
set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 15)

ylabel('Mean Firing Rate')
p_MFR = signrank(LADDER_MFR(:,1), LADDER_MFR(:,2))

% Rin LADDER

figure()
for u = 1:length(LADDER_Rin(:,1))
plot(x, LADDER_Rin(u,:),'-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'k')
hold on
end
means = [nanmean(LADDER_Rin(:,1)), nanmean(LADDER_Rin(:,2)), nanmean(LADDER_Rin(:,3))];
plot(x, means, '-md', 'MarkerSize', 20, 'MarkerFaceColor', 'm',  'LineWidth', 3)

xticks([1,2,3])
set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 15)
ylabel('Rin')
p_Rin = signrank(LADDER_Rin(:,1), LADDER_Rin(:,2))

% Vrm LADDER

figure()
for u = 1:length(LADDER_Vrm(:,1))
plot(x, LADDER_Vrm(u,:),'-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'k')
hold on
end
means = [nanmean(LADDER_Vrm(:,1)), nanmean(LADDER_Vrm(:,2)), nanmean(LADDER_Vrm(:,3))];
plot(x, means, '-md', 'MarkerSize', 20, 'MarkerFaceColor', 'm',  'LineWidth', 3)

xticks([1,2,3])
set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 15)

ylabel('Resting Membrane Potential')
p_Vrm = signrank(LADDER_Vrm(:,1), LADDER_Vrm(:,2))

% ADP LADDER

figure()
for u = 1:length(LADDER_ADP(:,1))
plot(x, LADDER_ADP(u,:),'-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'k')
hold on
end
means = [nanmean(LADDER_ADP(:,1)), nanmean(LADDER_ADP(:,2)), nanmean(LADDER_ADP(:,3))];
plot(x, means, '-md', 'MarkerSize', 20, 'MarkerFaceColor', 'm',  'LineWidth', 3)

xticks([1,2,3])
set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 15)

ylabel('ADP')
p_ADP = signrank(LADDER_ADP(:,1), LADDER_ADP(:,2))

% AHP LADDER

figure()
for u = 1:length(LADDER_AHP(:,1))
plot(x, LADDER_AHP(u,:),'-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'k')
hold on
end
means = [nanmean(LADDER_AHP(:,1)), nanmean(LADDER_AHP(:,2)), nanmean(LADDER_AHP(:,3))];
plot(x, means, '-md', 'MarkerSize', 20, 'MarkerFaceColor', 'm',  'LineWidth', 3)

xticks([1,2,3])
set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 15)

ylabel('AHP')
p_AHP = signrank(LADDER_AHP(:,1), LADDER_AHP(:,2))
end

LADDERS = struct;
LADDERS.AHP = LADDER_AHP;
LADDERS.ADP = LADDER_ADP;
LADDERS.Vrm  = LADDER_Vrm;
LADDERS.Rin = LADDER_Rin;
LADDERS.MFR = LADDER_MFR;
LADDERS.IFR = LADDER_IFR_ave;
LADDERS.count = count;
LADDERS.DATA = DATA;
end

