% In vivo plots, compare conditions: Hm4di + CNO vs CNO only
% Bottorff manuscript - updated Dec. 2022

clearvars -except MDCELLS CNO_MASTER *CELL* *LFP* statetimes

% CELL_DATA.deprived_ladder = MD_ladder_deprived_clean;
% CELL_DATA.control_ladder = MD_ladder_control_clean;
% CELL_DATA.ladder_period = period_v2;
% CELL_DATA.ladder_MD2end = MD2_end_v2;
% CELL_DATA.ladder_MD4end = MD4_end_v2;
% CELL_DATA.bintime = bintime2;
% CELL_DATA.FR_deprived = all_FR_normalized_RSU;
% CELL_DATA.FR_control = FR_normalized_RSU_control;

Hm4di_CNO_MDCELLS = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/TlabDocs/ACh_Manuscript/PaperFigs_drafts_2022/JB_MASTER_Hm4di_CNO_CELLS_v5.mat');
CNOonly_MDCELLS = load('/Users/Juliet/Library/CloudStorage/GoogleDrive-jbottorff@brandeis.edu/My Drive/TlabDocs/ACh_Manuscript/PaperFigs_drafts_2022/JB_MASTER_CNOonly_CELLS_v2.mat');

CELL_DATA_Hm4di = inVivoAnalysis_BFACh_JBmanuscript(Hm4di_CNO_MDCELLS.MDCELLS, 0);
CELL_DATA_CNOonly = inVivoAnalysis_BFACh_JBmanuscript(CNOonly_MDCELLS.MDCELLS, 0);


% plot box plot of percent change from baseline for deprived and control,
% BFinhib vs CNO only
Colors = [0,0,0;1,0,1; 0.5,0.5,0.5; 0,0,1];

BFinhib_deprivedDiff = CELL_DATA_Hm4di.deprived_ladder(:,3)./CELL_DATA_Hm4di.deprived_ladder(:,1);
BFinhib_controlDiff = CELL_DATA_Hm4di.control_ladder(:,3)./CELL_DATA_Hm4di.control_ladder(:,1);

CNOonly_deprivedDiff = CELL_DATA_CNOonly.deprived_ladder(:,3)./CELL_DATA_CNOonly.deprived_ladder(:,1);
CNOonly_controlDiff = CELL_DATA_CNOonly.control_ladder(:,3)./CELL_DATA_CNOonly.control_ladder(:,1);

figure()
data_toplot_1 = padmat((CNOonly_controlDiff.*100), (CNOonly_deprivedDiff.*100),2);
data_toplot_2 = padmat((BFinhib_controlDiff.*100), (BFinhib_deprivedDiff.*100),2);
data_toplot = padmat(data_toplot_1, data_toplot_2, 2);
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);
hold on
plot([0,4.5], [100,100], 'k--')
box off
ylabel('Percent of Baseline FR')
set(gca,'XTickLabel',{'Cont. Hem CNOonly', 'Dep. Hem CNOonly', 'Cont. Hem Hm4di+CNO', 'Dep. Hem Hm4di+CNO'},'XTickLabelRotation',45,'FontSize',15);
 
%ylim([0 300])
xticks([1 2 3 4])
p_ranksum_PercChange = ranksum(data_toplot(:,1), data_toplot(:,2))
[p_kruskalwallis_percChange, tbl, stats] = kruskalwallis(data_toplot)
multcompare(stats)

Colors_Face = [1,1,1; 1,0,1; 1,1,1; 0,0,1];
Colors = [1,0,1;1,0,1; 0,0,1; 0,0,1];
% plot MD2 vs MD4 for each deprived group
BFinhib_MD2Diff = CELL_DATA_Hm4di.deprived_ladder(:,2)./CELL_DATA_Hm4di.deprived_ladder(:,1);
BFinhib_MD4Diff = CELL_DATA_Hm4di.deprived_ladder(:,3)./CELL_DATA_Hm4di.deprived_ladder(:,1);

CNOonly_MD2Diff = CELL_DATA_CNOonly.deprived_ladder(:,2)./CELL_DATA_CNOonly.deprived_ladder(:,1);
CNOonly_MD4Diff = CELL_DATA_CNOonly.deprived_ladder(:,3)./CELL_DATA_CNOonly.deprived_ladder(:,1);

figure()
data_toplot_1 = padmat((CNOonly_MD2Diff.*100), (CNOonly_MD4Diff.*100),2);
data_toplot_2 = padmat((BFinhib_MD2Diff.*100), (BFinhib_MD4Diff.*100),2);
data_toplot = padmat(data_toplot_1, data_toplot_2, 2);
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors_Face,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);
hold on
plot([0,4.5], [100,100], 'k--')
box off
ylabel('Percent of Baseline FR')
set(gca,'XTickLabel',{'CNOonly MD2', 'CNOonly MD4', 'Hm4di+CNO MD2', 'Hm4di+CNO MD4'},'XTickLabelRotation',45,'FontSize',15);
title ('Deprived Hem')
%ylim([0 300])
xticks([1 2 3 4])
p_ranksum_PercChange = ranksum(data_toplot(:,1), data_toplot(:,2))
[p_kruskalwallis_percChange, tbl, stats] = kruskalwallis(data_toplot)
multcompare(stats)


% plot baseline vs MD2 and baseline vs MD4 FR's for hm4di+CNO and CNO only

figure()
hold on
for i = 1:length(CELL_DATA_Hm4di.deprived_ladder(:,1))
    plot(CELL_DATA_Hm4di.deprived_ladder(i,1), CELL_DATA_Hm4di.deprived_ladder(i,2), 'bo', 'MarkerFaceColor', [1,1,1], 'LineWidth', 3)
end

for i = 1:length(CELL_DATA_CNOonly.deprived_ladder(:,1))
    plot(CELL_DATA_CNOonly.deprived_ladder(i,1), CELL_DATA_CNOonly.deprived_ladder(i,2), 'mo', 'MarkerFaceColor', [1,1,1], 'LineWidth', 3)
end
plot([10e-4, 10e2], [10e-4, 10e2], 'k', 'LineWidth', 3)
xlim([10e-4, 10e2])
ylim([10e-4, 10e2])
box off
set(gca,'FontSize',15);
set(gca, 'YScale', 'log', 'XScale', 'log')
xlabel('Baseline Firing Rate (Hz)')
ylabel('MD2 Firing Rate (Hz)')

figure()
hold on
for i = 1:length(CELL_DATA_Hm4di.deprived_ladder(:,1))
    plot(CELL_DATA_Hm4di.deprived_ladder(i,1), CELL_DATA_Hm4di.deprived_ladder(i,3), 'bo', 'MarkerFaceColor', 'b', 'LineWidth', 3)
end
for i = 1:length(CELL_DATA_CNOonly.deprived_ladder(:,1))
    plot(CELL_DATA_CNOonly.deprived_ladder(i,1), CELL_DATA_CNOonly.deprived_ladder(i,3), 'mo', 'MarkerFaceColor', 'm', 'LineWidth', 3)
end
plot([10e-4, 10e2], [10e-4, 10e2], 'k', 'LineWidth', 3)
xlim([10e-4, 10e2])
ylim([10e-4, 10e2])
box off
set(gca,'FontSize',15);
set(gca, 'YScale', 'log', 'XScale', 'log')
xlabel('Baseline Firing Rate (Hz)')
ylabel('MD4 Firing Rate (Hz)')


CNOonly_meanFR_deprived = nanmean(CELL_DATA_CNOonly.FR_deprived);
Hm4di_meanFR_deprived = nanmean(CELL_DATA_Hm4di.FR_deprived);

CNOonly_semFR_deprived=(nanstd(CELL_DATA_CNOonly.FR_deprived))./(sqrt(length(CELL_DATA_CNOonly.FR_deprived(:,1)))-1);
Hm4di_semFR_deprived=(nanstd(CELL_DATA_Hm4di.FR_deprived))./(sqrt(length(CELL_DATA_Hm4di.FR_deprived(:,1)))-1);

if length(CNOonly_meanFR_deprived) < length(Hm4di_meanFR_deprived)
    Hm4di_meanFR_deprived = Hm4di_meanFR_deprived(1:length(CNOonly_meanFR_deprived));
    Hm4di_semFR_deprived = Hm4di_semFR_deprived(1:length(CNOonly_semFR_deprived));
else
    CNOonly_meanFR_deprived = CNOonly_meanFR_deprived(1:length(Hm4di_meanFR_deprived));
    CNOonly_semFR_deprived = CNOonly_semFR_deprived(1:length(Hm4di_semFR_deprived));
end

xaxis = 1:length(Hm4di_meanFR_deprived);
