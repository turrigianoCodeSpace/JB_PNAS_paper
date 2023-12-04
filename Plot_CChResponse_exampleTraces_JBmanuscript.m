% Plot example traces for CCh response experiments
% ACh manuscript Figure 4E

DATA_ctrl = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/ALL_ctrl_CCh_Analysis.mat');
DATA_m1_3 = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/m1_3sh_CCh_Analysis.mat');


DATA_ctrl = DATA_ctrl.FI_DATA;
DATA_m1_3 = DATA_m1_3.FI_DATA;

Colors = [0 0 0; 0 1 1; 0 0 1];

%CHOOSE which cells from each group for example: 
ctrl = 5; 
m1_3 = 16; 

baseline_trace = 10;
CCh_trace = 136;
wash_trace = 248;

current_step = 250;

figure(1)
plot(DATA_ctrl.aDAT{ctrl}(:,baseline_trace), 'Color', Colors(1,:), 'LineWidth', 2.5)
hold on
plot([12000, 14000], [-30,-30], 'k', 'LineWidth', 2.5) % scale bar: 0.2 sec
plot([12000, 12000], [-30, 10], 'k', 'LineWidth', 2.5) % scale bar: 40 mV
plot([0, 30000], [mean(DATA_ctrl.aDAT{ctrl}(2000:4000,baseline_trace)), mean(DATA_ctrl.aDAT{ctrl}(2000:4000,baseline_trace))], 'k--')
set(gca, 'box', 'off', 'visible', 'off')
ylim([-80, 60])
xlim([3000, 20000])
title('Ctrl baseline')

figure(11)
plot(DATA_ctrl.aDAT{ctrl}(:,baseline_trace), 'Color', Colors(1,:), 'LineWidth', 2.5)
hold on
plot([12000, 14000], [-72,-72], 'k', 'LineWidth', 2.5) % scale bar: 0.2 sec
plot([12000, 12000], [-72, -71], 'k', 'LineWidth', 2.5) % scale bar: 1 mV
plot([0, 30000], [mean(DATA_ctrl.aDAT{ctrl}(2000:4000,baseline_trace)), mean(DATA_ctrl.aDAT{ctrl}(2000:4000,baseline_trace))], 'k--')
set(gca, 'box', 'off', 'visible', 'off')
ylim([-75, -71])
xlim([10000, 15000])
title('Ctrl baseline')

figure(2)
plot(DATA_ctrl.aDAT{ctrl}(:,CCh_trace), 'Color', Colors(1,:), 'LineWidth', 2.5)
hold on
plot([0, 30000], [mean(DATA_ctrl.aDAT{ctrl}(2000:4000,CCh_trace)), mean(DATA_ctrl.aDAT{ctrl}(2000:4000,CCh_trace))], 'k--')
set(gca, 'box', 'off', 'visible', 'off')
ylim([-80, 60])
xlim([3000, 20000])
title('Ctrl CCh')

figure(12)
plot(DATA_ctrl.aDAT{ctrl}(:,CCh_trace), 'Color', Colors(1,:), 'LineWidth', 2.5)
hold on
plot([12000, 14000], [-67,-67], 'k', 'LineWidth', 2.5) % scale bar: 0.2 sec
plot([12000, 12000], [-67, -66], 'k', 'LineWidth', 2.5) % scale bar: 1 mV
plot([0, 30000], [mean(DATA_ctrl.aDAT{ctrl}(2000:4000,CCh_trace)), mean(DATA_ctrl.aDAT{ctrl}(2000:4000,CCh_trace))], 'k--')
set(gca, 'box', 'off', 'visible', 'off')
ylim([-71, -66])
xlim([10000, 15000])
title('Ctrl CCh')

figure(3)
plot(DATA_ctrl.aDAT{ctrl}(:,wash_trace), 'Color', Colors(1,:), 'LineWidth', 2.5)
hold on
plot([0, 30000], [mean(DATA_ctrl.aDAT{ctrl}(2000:4000,wash_trace)), mean(DATA_ctrl.aDAT{ctrl}(2000:4000,wash_trace))], 'k--')
set(gca, 'box', 'off', 'visible', 'off')
ylim([-80, 60])
xlim([3000, 20000])
title('Ctrl wash')

figure(13)
plot(DATA_ctrl.aDAT{ctrl}(:,wash_trace), 'Color', Colors(1,:), 'LineWidth', 2.5)
hold on
plot([12000, 14000], [-72,-72], 'k', 'LineWidth', 2.5) % scale bar: 0.2 sec
plot([12000, 12000], [-72, -71], 'k', 'LineWidth', 2.5) % scale bar: 1 mV
plot([0, 30000], [mean(DATA_ctrl.aDAT{ctrl}(2000:4000,wash_trace)), mean(DATA_ctrl.aDAT{ctrl}(2000:4000,wash_trace))], 'k--')
set(gca, 'box', 'off', 'visible', 'off')
ylim([-75, -71])
xlim([10000, 15000])
title('Ctrl wash')

figure(7)
plot(DATA_m1_3.aDAT{m1_3}(:,baseline_trace), 'Color', Colors(3,:), 'LineWidth', 2.5)
hold on
plot([0, 30000], [mean(DATA_m1_3.aDAT{m1_3}(2000:4000,baseline_trace)), mean(DATA_m1_3.aDAT{m1_3}(2000:4000,baseline_trace))], 'k--')
set(gca, 'box', 'off', 'visible', 'off')
xlim([3000, 20000])
ylim([-80, 60])
title('m1_3 baseline')

figure(17)
plot(DATA_m1_3.aDAT{m1_3}(:,baseline_trace), 'Color', Colors(3,:), 'LineWidth', 2.5)
hold on
plot([12000, 14000], [-74,-74], 'k', 'LineWidth', 2.5) % scale bar: 0.2 sec
plot([12000, 12000], [-74, -73], 'k', 'LineWidth', 2.5) % scale bar: 1 mV
plot([0, 30000], [mean(DATA_m1_3.aDAT{m1_3}(2000:4000,baseline_trace)), mean(DATA_m1_3.aDAT{m1_3}(2000:4000,baseline_trace))], 'k--')
set(gca, 'box', 'off', 'visible', 'off')
ylim([-76.7, -72.7])
xlim([10000, 15000])
title('m1_3 baseline')

figure(8)
plot(DATA_m1_3.aDAT{m1_3}(:,CCh_trace), 'Color', Colors(3,:), 'LineWidth', 2.5)
hold on
plot([0, 30000], [mean(DATA_m1_3.aDAT{m1_3}(2000:4000,CCh_trace)), mean(DATA_m1_3.aDAT{m1_3}(2000:4000,CCh_trace))], 'k--')
set(gca, 'box', 'off', 'visible', 'off')
xlim([3000, 20000])
ylim([-80, 60])
title('m1_3 CCh')

figure(18)
plot(DATA_m1_3.aDAT{m1_3}(:,CCh_trace), 'Color', Colors(3,:), 'LineWidth', 2.5)
hold on
plot([12000, 14000], [-74,-74], 'k', 'LineWidth', 2.5) % scale bar: 0.2 sec
plot([12000, 12000], [-74, -73], 'k', 'LineWidth', 2.5) % scale bar: 1 mV
plot([0, 30000], [mean(DATA_m1_3.aDAT{m1_3}(2000:4000,CCh_trace)), mean(DATA_m1_3.aDAT{m1_3}(2000:4000,CCh_trace))], 'k--')
set(gca, 'box', 'off', 'visible', 'off')
ylim([-77, -73])
xlim([10000, 15000])
title('m1_3 CCh')

figure(9)
plot(DATA_m1_3.aDAT{m1_3}(:,wash_trace), 'Color', Colors(3,:), 'LineWidth', 2.5)
hold on
plot([0, 30000], [mean(DATA_m1_3.aDAT{m1_3}(2000:4000,wash_trace)), mean(DATA_m1_3.aDAT{m1_3}(2000:4000,wash_trace))], 'k--')
set(gca, 'box', 'off', 'visible', 'off')
xlim([3000, 20000])
ylim([-80, 60])
title('m1_3 wash')

figure(19)
plot(DATA_m1_3.aDAT{m1_3}(:,wash_trace), 'Color', Colors(3,:), 'LineWidth', 2.5)
hold on
plot([12000, 14000], [-72,-72], 'k', 'LineWidth', 2.5) % scale bar: 0.2 sec
plot([12000, 12000], [-72, -71], 'k', 'LineWidth', 2.5) % scale bar: 1 mV
plot([0, 30000], [mean(DATA_m1_3.aDAT{m1_3}(2000:4000,wash_trace)), mean(DATA_m1_3.aDAT{m1_3}(2000:4000,wash_trace))], 'k--')
set(gca, 'box', 'off', 'visible', 'off')
ylim([-75, -71])
xlim([10000, 15000])
title('m1_3 wash')



figure(10)
plot([0, 5000], [0.5, 0.5], 'k', 'LineWidth', 2.5)
hold on
plot([5000, 5000], [0.5, 1], 'k', 'LineWidth', 2.5)
plot([5000, 15000], [1, 1], 'k', 'LineWidth', 2.5)
plot([15000, 15000], [0.5, 1], 'k', 'LineWidth', 2.5)
plot([15000, 20000], [0.5, 0.5], 'k', 'LineWidth', 2.5)
set(gca, 'box', 'off', 'visible', 'off')
ylim([0, 4])
