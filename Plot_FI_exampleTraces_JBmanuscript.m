% Plot example traces for FI curves

shRNA = 0; % are you plotting shRNA experiments? 1 = Figure 5; 0 = Figure 3

if ~shRNA
FI_Hm4di = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/101921thru110521_Hm4di_FI_DATA_v2.mat'); % Hm4di in V1m only, 24 hours CNO!
FI_CNOonly = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/120222thru121922_EF1a_emptyVector_FI_DATA.mat'); % EF1a empty vector in V1m + 24 hours CNO!
FI_BF_Hm4di = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/111821thru112421_Hm4di_FI_DATA_v2.mat'); % Hm4di in V1m and BF, 24 hours CNO!
FI_BF_CNOonly = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/111821thru112421_CNOonly_FI_DATA_v2.mat'); % Hm4di in V1m and BF, 24 hours CNO!
else
%FI_BF_Hm4di = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/090822thru101822_m1_3shRNA_Hm4di_FI_DATA.mat'); % Hm4di + m1-3shRNA in V1m combined + CNO, 14+ days knockdown
%FI_BF_CNOonly = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/090822thru121922_m1_3shRNA_CNOonly_FI_DATA.mat'); % m1-3shRNA in V1m combined + CNO, 14+ days knockdown
FI_BF_Hm4di = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/072423thru072823_m1sh_Hm4di_FI_DATA.mat'); % m1-3shRNA and Hm4di in V1m and 24 hours CNO - 6-10 days knockdown (p24-p28)
FI_BF_CNOonly = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/Mini_Data/072423thru072823_m1sh_FI_DATA.mat'); % m1-3shRNA in V1m and 24 hours CNO - 6-10 days knockdown (p24-p28)
end

if ~shRNA
FI_Hm4di = FI_Hm4di.FI_DATA;
FI_CNOonly = FI_CNOonly.FI_DATA;
end
FI_BF_Hm4di = FI_BF_Hm4di.FI_DATA;
FI_BF_CNOonly = FI_BF_CNOonly.FI_DATA;

Colors = [0 0 0; 1 0 1; 0.5 0.5 0.5; 0.58 0.82 0.98];

%CHOOSE which cells from each group for example: 
if ~shRNA
Hm4di_cell = 10;
CNOonly_cell = 25;
BF_Hm4di_cell = 16; 
BF_CNOonly_cell = 5; 

else
% SHORTER KD (6-10D) GROUPS!
BF_CNOonly_cell = 7;
BF_Hm4di_cell = 11;
% LONGER KD GROUPS (14+ DAYS)
%BF_Hm4di_cell = 5; % for shRNA expts; 
%BF_CNOonly_cell = 29; % for shRNA expts; 
end

current_step = 200;
num_curr_inj = current_step/20; 

if ~shRNA
figure(1)
plot(FI_Hm4di.aDAT{Hm4di_cell}(:,num_curr_inj), 'Color', Colors(2,:), 'LineWidth', 2.5)
set(gca, 'box', 'off', 'visible', 'off')
ylim([-80, 60])
xlim([5000, 25000])
title('Hm4di')

figure(2)
plot(FI_CNOonly.aDAT{CNOonly_cell}(:,num_curr_inj), 'Color', Colors(1,:), 'LineWidth', 2.5)
hold on
plot([6000, 8000], [-30,-30], 'k', 'LineWidth', 2.5) % scale bar: 0.2 sec
plot([6000, 6000], [-30, 10], 'k', 'LineWidth', 2.5) % scale bar: 40 mV
set(gca, 'box', 'off', 'visible', 'off')
xlim([5000, 25000])
ylim([-80, 60])
title('CNO only')
end

figure(3)
plot(FI_BF_Hm4di.aDAT{BF_Hm4di_cell}(:,num_curr_inj), 'Color', Colors(4,:), 'LineWidth', 2.5)
set(gca, 'box', 'off', 'visible', 'off')
xlim([5000, 25000])
ylim([-80, 60])
title('BF Hm4di')

figure(4)
plot(FI_BF_CNOonly.aDAT{BF_CNOonly_cell}(:,num_curr_inj), 'Color', Colors(3,:), 'LineWidth', 2.5)
set(gca, 'box', 'off', 'visible', 'off')
xlim([5000, 25000])
ylim([-80, 60])
title('BF CNO only')

figure(5)
plot([0, 5000], [0.5, 0.5], 'k', 'LineWidth', 2.5)
hold on
plot([5000, 5000], [0.5, 1], 'k', 'LineWidth', 2.5)
plot([5000, 15000], [1, 1], 'k', 'LineWidth', 2.5)
plot([15000, 15000], [0.5, 1], 'k', 'LineWidth', 2.5)
plot([15000, 20000], [0.5, 0.5], 'k', 'LineWidth', 2.5)
set(gca, 'box', 'off', 'visible', 'off')
ylim([0, 4])
