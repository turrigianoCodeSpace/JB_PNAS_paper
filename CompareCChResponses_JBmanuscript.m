% Compare Response to CCh in ctrl vs M1 knockdown cells
% Calls Script "CChResponse_Plots_JBmanuscript.m"
% Uses ladder plot data from that script to compare changes from baseline
% to CCh application
% ACh manuscript Figure 4F

clearvars


LADDERS_ctrl = CChResponse_Plots_JBmanuscript('ctrl', 0);
LADDERS_m1_3 = CChResponse_Plots_JBmanuscript('m1_3', 0);

close all

Colors  = [0,0,0; 0,0,1; 0,1,1];

% available ladders
%LADDERS.AHP = LADDER_AHP;
%LADDERS.ADP = LADDER_ADP;
%LADDERS.Vrm  = LADDER_Vrm;
%LADDERS.Rin = LADDER_Rin;
%LADDERS.MFR = LADDER_MFR;
%LADDERS.IFR = LADDER_IFR_ave;
%LADDERS.count = count;
%LADDERS.DATA = DATA;


AHP_diff_ctrl = LADDERS_ctrl.AHP(:,2) - LADDERS_ctrl.AHP(:,1);
ADP_diff_ctrl = LADDERS_ctrl.ADP(:,2) - LADDERS_ctrl.ADP(:,1);
Vrm_diff_ctrl = LADDERS_ctrl.Vrm(:,2) - LADDERS_ctrl.Vrm(:,1);
Rin_diff_ctrl = LADDERS_ctrl.Rin(:,2) - LADDERS_ctrl.Rin(:,1);
MFR_diff_ctrl = LADDERS_ctrl.MFR(:,2) - LADDERS_ctrl.MFR(:,1);
IFR_diff_ctrl = LADDERS_ctrl.IFR(:,2) - LADDERS_ctrl.IFR(:,1);

AHP_diff_m1_3 = LADDERS_m1_3.AHP(:,2) - LADDERS_m1_3.AHP(:,1);
ADP_diff_m1_3 = LADDERS_m1_3.ADP(:,2) - LADDERS_m1_3.ADP(:,1);
Vrm_diff_m1_3 = LADDERS_m1_3.Vrm(:,2) - LADDERS_m1_3.Vrm(:,1);
Rin_diff_m1_3 = LADDERS_m1_3.Rin(:,2) - LADDERS_m1_3.Rin(:,1);
MFR_diff_m1_3 = LADDERS_m1_3.MFR(:,2) - LADDERS_m1_3.MFR(:,1);
IFR_diff_m1_3 = LADDERS_m1_3.IFR(:,2) - LADDERS_m1_3.IFR(:,1);


% Plot AHP differences
figure()

data_toplot = padmat(AHP_diff_ctrl, AHP_diff_m1_3, 2);
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel('AHP CCh - Baseline')

set(gca,'XTickLabel',{'Ctrl', 'M1-3sh'},'XTickLabelRotation',45,'FontSize',15);
xticks([1 2])
%ylim([0 30])
p_ranksum_AHP_m1_3 = ranksum(data_toplot(:,1), data_toplot(:,2))



% Plot ADP differences
figure()

data_toplot = padmat(ADP_diff_ctrl, ADP_diff_m1_3, 2);
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel('ADP CCh - Baseline')

set(gca,'XTickLabel',{'Ctrl', 'M1-3sh'},'XTickLabelRotation',45,'FontSize',15);
xticks([1 2])
%ylim([0 30])
p_ranksum_ADP_m1_3 = ranksum(data_toplot(:,1), data_toplot(:,2))


% Plot Rin differences
figure()

data_toplot = padmat(Rin_diff_ctrl, Rin_diff_m1_3, 2);
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel('Rin CCh - Baseline')

set(gca,'XTickLabel',{'Ctrl', 'M1-3sh'},'XTickLabelRotation',45,'FontSize',15);
xticks([1 2])
p_ranksum_Rin_m1_3 = ranksum(data_toplot(:,1), data_toplot(:,2))

% Plot MFR differences
figure()

data_toplot = padmat(MFR_diff_ctrl, MFR_diff_m1_3, 2);
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel('MFR CCh - Baseline')

set(gca,'XTickLabel',{'Ctrl', 'M1-3sh'},'XTickLabelRotation',45,'FontSize',15);
xticks([1 2])
p_ranksum_MFR_m1_3 = ranksum(data_toplot(:,1), data_toplot(:,2))


% Plot IFR differences
figure()

data_toplot = padmat(IFR_diff_ctrl, IFR_diff_m1_3, 2);
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel('IFR CCh - Baseline')

set(gca,'XTickLabel',{'Ctrl', 'M1-3sh'},'XTickLabelRotation',45,'FontSize',15);
xticks([1 2])
p_ranksum_IFR_m1_3 = ranksum(data_toplot(:,1), data_toplot(:,2))


% LADDER COMPARISONS for ctrl vs m1-3
Colors  = [0,0,0; 0,1,1; 0,0,1];

names = {'ctrl baseline', 'ctrl CCh', 'ctrl wash', 'm1-3 baseline', 'm1-3 CCh', 'm1-3 wash'};

% AHP LADDER

figure()
delta_LADDERS_ctrl.AHP = LADDERS_ctrl.AHP-LADDERS_ctrl.AHP(:,1);
delta_LADDERS_m1_3.AHP = LADDERS_m1_3.AHP-LADDERS_m1_3.AHP(:,1);
for u = 1:length(delta_LADDERS_ctrl.AHP(:,1))
plot([1,2,3], delta_LADDERS_ctrl.AHP(u,1:3),'-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'k')
hold on
end
means = [nanmean(delta_LADDERS_ctrl.AHP(:,1)), nanmean(delta_LADDERS_ctrl.AHP(:,2)), nanmean(delta_LADDERS_ctrl.AHP(:,3))];
plot([1,2,3], means, '-kd', 'MarkerSize', 20, 'MarkerFaceColor', [.5,.5,.5],  'LineWidth', 3)

for u = 1:length(delta_LADDERS_m1_3.AHP(:,1))
plot([4,5,6], delta_LADDERS_m1_3.AHP(u,1:3),'-bo', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'b')
hold on
end
means = [nanmean(delta_LADDERS_m1_3.AHP(:,1)), nanmean(delta_LADDERS_m1_3.AHP(:,2)), nanmean(delta_LADDERS_m1_3.AHP(:,3))];
plot([4,5,6], means, '-bd', 'MarkerSize', 20, 'MarkerFaceColor', 'b',  'LineWidth', 3)
xticks([1,2,3,4,5,6])
set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 15, 'XTickLabelRotation', 45)
ylabel('AHP')
xlim([0.5, 6.5])

% ADP LADDER
figure()
for u = 1:length(LADDERS_ctrl.ADP(:,1))
plot([1,2,3], LADDERS_ctrl.ADP(u,1:3),'-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'k')
hold on
end
means = [nanmean(LADDERS_ctrl.ADP(:,1)), nanmean(LADDERS_ctrl.ADP(:,2)), nanmean(LADDERS_ctrl.ADP(:,3))];
plot([1,2,3], means, '-kd', 'MarkerSize', 20, 'MarkerFaceColor', [.5,.5,.5],  'LineWidth', 3)

for u = 1:length(LADDERS_m1_3.ADP(:,1))
plot([4,5,6], LADDERS_m1_3.ADP(u,1:3),'-bo', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'b')
hold on
end
means = [nanmean(LADDERS_m1_3.ADP(:,1)), nanmean(LADDERS_m1_3.ADP(:,2)), nanmean(LADDERS_m1_3.ADP(:,3))];
plot([4,5,6], means, '-bd', 'MarkerSize', 20, 'MarkerFaceColor', 'b',  'LineWidth', 3)

xticks([1,2,3,4,5,6])
set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 15, 'XTickLabelRotation', 45)
ylabel('ADP')
xlim([0.5, 6.5])


% Rin LADDER
delta_LADDERS_ctrl.Rin = LADDERS_ctrl.Rin-LADDERS_ctrl.Rin(:,1);
delta_LADDERS_m1_3.Rin = LADDERS_m1_3.Rin-LADDERS_m1_3.Rin(:,1);
figure()
for u = 1:length(delta_LADDERS_ctrl.Rin(:,1))
plot([1,2,3], delta_LADDERS_ctrl.Rin(u,1:3),'-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'k')
hold on
end
means = [nanmean(delta_LADDERS_ctrl.Rin(:,1)), nanmean(delta_LADDERS_ctrl.Rin(:,2)), nanmean(delta_LADDERS_ctrl.Rin(:,3))];
plot([1,2,3], means, '-kd', 'MarkerSize', 20, 'MarkerFaceColor', [.5,.5,.5],  'LineWidth', 3)

for u = 1:length(delta_LADDERS_m1_3.Rin(:,1))
plot([4,5,6], delta_LADDERS_m1_3.Rin(u,1:3),'-bo', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'b')
hold on
end
means = [nanmean(delta_LADDERS_m1_3.Rin(:,1)), nanmean(delta_LADDERS_m1_3.Rin(:,2)), nanmean(delta_LADDERS_m1_3.Rin(:,3))];
plot([4,5,6], means, '-bd', 'MarkerSize', 20, 'MarkerFaceColor', 'b',  'LineWidth', 3)

xticks([1,2,3,4,5,6])
set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 15, 'XTickLabelRotation', 45)
ylabel('Rin')
xlim([0.5, 6.5])

% MFR LADDER
delta_LADDERS_ctrl.MFR = LADDERS_ctrl.MFR-LADDERS_ctrl.MFR(:,1);
delta_LADDERS_m1_3.MFR = LADDERS_m1_3.MFR-LADDERS_m1_3.MFR(:,1);
figure()
for u = 1:length(delta_LADDERS_ctrl.MFR(:,1))
plot([1,2,3], delta_LADDERS_ctrl.MFR(u,1:3),'-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'k')
hold on
end
means = [nanmean(delta_LADDERS_ctrl.MFR(:,1)), nanmean(delta_LADDERS_ctrl.MFR(:,2)), nanmean(delta_LADDERS_ctrl.MFR(:,3))];
plot([1,2,3], means, '-kd', 'MarkerSize', 20, 'MarkerFaceColor', [.5,.5,.5],  'LineWidth', 3)

for u = 1:length(delta_LADDERS_m1_3.MFR(:,1))
plot([4,5,6], delta_LADDERS_m1_3.MFR(u,1:3),'-bo', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'b')
hold on
end
means = [nanmean(delta_LADDERS_m1_3.MFR(:,1)), nanmean(delta_LADDERS_m1_3.MFR(:,2)), nanmean(delta_LADDERS_m1_3.MFR(:,3))];
plot([4,5,6], means, '-bd', 'MarkerSize', 20, 'MarkerFaceColor', 'b',  'LineWidth', 3)

xticks([1,2,3,4,5,6])
set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 15, 'XTickLabelRotation', 45)
ylabel('MFR')
xlim([0.5, 6.5])

% IFR LADDER
delta_LADDERS_ctrl.IFR = LADDERS_ctrl.IFR-LADDERS_ctrl.IFR(:,1);
delta_LADDERS_m1_3.IFR = LADDERS_m1_3.IFR-LADDERS_m1_3.IFR(:,1);
figure()
for u = 1:length(delta_LADDERS_ctrl.IFR(:,1))
plot([1,2,3], delta_LADDERS_ctrl.IFR(u,1:3),'-ko', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'k')
hold on
end
means = [nanmean(delta_LADDERS_ctrl.IFR(:,1)), nanmean(delta_LADDERS_ctrl.IFR(:,2)), nanmean(delta_LADDERS_ctrl.IFR(:,3))];
plot([1,2,3], means, '-kd', 'MarkerSize', 20, 'MarkerFaceColor', [.5,.5,.5],  'LineWidth', 3)

for u = 1:length(delta_LADDERS_m1_3.IFR(:,1))
plot([4,5,6], delta_LADDERS_m1_3.IFR(u,1:3),'-bo', 'MarkerSize', 13, 'MarkerFaceColor', 'w', 'LineWidth', 3, 'MarkerEdgeColor', 'b')
hold on
end
means = [nanmean(delta_LADDERS_m1_3.IFR(:,1)), nanmean(delta_LADDERS_m1_3.IFR(:,2)), nanmean(delta_LADDERS_m1_3.IFR(:,3))];
plot([4,5,6], means, '-bd', 'MarkerSize', 20, 'MarkerFaceColor', 'b',  'LineWidth', 3)

xticks([1,2,3,4,5,6])
set (gca, 'xticklabel', names, 'box', 'off', 'fontsize', 15, 'XTickLabelRotation', 45)
ylabel('IFR')
xlim([0.5, 6.5])
