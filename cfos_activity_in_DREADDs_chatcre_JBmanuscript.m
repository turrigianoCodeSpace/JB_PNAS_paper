% plot cfos data in DREADDs chatcre animals from Sydney Padgett
% ACh manuscript Figure 1D


% Animal	CNO Condition	Sleep Dep Time	% of RFP cells cfos+	
% SP1	No CNO	7:30 AM	18.92%	
% SP5	No CNO	7:30 AM	21.74%	
% SP53  No CNO  7:30 AM 30.3%
% SP54  No CNO  7:30 AM 15.15%
% SP2	 48 hrs 	7:30 AM	10.194%	
% SP3	 48 hrs 	7:30 AM	13.619%	
% SP4	 48 hrs 	7:30 AM	7.965%	
% SP55   48 hrs     7:30 AM 12.903%
% SP66   48 hrs     7:30 AM 4.585%
% SP7	 12 hrs  	7:30 AM	9.851%	
% SP8	 12 hrs  	7:30 AM	12.821%	
% SP11	 12 hrs  	7:30 AM	5.597%	
% SP65   12 hrs     7:30 AM 13.818%
% SP67   12 hrs     7:30 AM 8.252%

% animals with no DREADDs: just ChAT and cFos
% SP27	7:30 AM	19.565%
% SP28	7:30 AM	15.948%
% SP29	7:30 AM	22.798%
% SP30	7:30 AM	16.816%
% SP33	7:30 AM	24.116%
% SP34	7:30 AM	16.389%


AM_NOCNO = [18.92, 21.74, 30.3, 15.15];
AM_48hrsCNO = [10.194, 13.619, 7.965, 12.903, 4.585];
AM_12hrsCNO = [9.851, 12.821, 5.597, 13.818, 8.252];

AM_NOCNO_NODREADDs = [19.565, 15.948, 22.798, 16.816, 24.116, 16.389];

mean_AMConditions = [nanmean([AM_NOCNO_NODREADDs, AM_NOCNO]), nanmean(AM_12hrsCNO), nanmean(AM_48hrsCNO)];
sem_AMConditions = [nanstd([AM_NOCNO_NODREADDs, AM_NOCNO])./sqrt(length([AM_NOCNO_NODREADDs, AM_NOCNO])-1), nanstd(AM_12hrsCNO)./sqrt(length(AM_12hrsCNO)-1), nanstd(AM_48hrsCNO)./sqrt(length(AM_48hrsCNO)-1)];

deviation = 0.35;
NOCNO_jitter = 1 + rand(1,length([AM_NOCNO_NODREADDs, AM_NOCNO])).*deviation - deviation/2;
AM_12hrs_jitter = 2 + rand(1,length(AM_12hrsCNO)).*deviation - deviation/2;
AM_48hrs_jitter = 3 + rand(1,length(AM_48hrsCNO)).*deviation - deviation/2;

figure()
b = bar(mean_AMConditions, 'EdgeColor', 'k', 'LineWidth', 5);
hold on
errorbar(mean_AMConditions, sem_AMConditions, 'k+')
plot(NOCNO_jitter(1:length(AM_NOCNO_NODREADDs)), AM_NOCNO_NODREADDs, 'x', 'Color', [0.3, 0.3, 0.3], 'LineWidth', 5, 'MarkerSize', 10)
plot(NOCNO_jitter(length(AM_NOCNO_NODREADDs)+1:end), AM_NOCNO, 'ko', 'LineWidth', 5, 'MarkerSize', 10)
plot(AM_12hrs_jitter, AM_12hrsCNO, 'ko', 'LineWidth', 5, 'MarkerSize', 10)
plot(AM_48hrs_jitter, AM_48hrsCNO, 'ko', 'LineWidth', 5, 'MarkerSize', 10)
xticks([1:3])
xlim([0.5,3.5])
set(gca, 'box', 'off', 'xticklabel', {'NO CNO', '12 hours CNO', '48 hours CNO'}, 'FontSize', 15)
b.FaceColor = 'flat';
b.CData(1,:) = [.7 .7 .7];
b.CData(2,:) = [.8 0 .1];
b.CData(3,:) = [.4 0 .1];
xtickangle(45)
ylabel('cFos+RFP+ cells / total RFP+ cells')
title ('Effect of Hm4di+CNO on ACh neuron activity')
DATA = padmat([AM_NOCNO_NODREADDs, AM_NOCNO]', AM_12hrsCNO',2);
DATA = padmat(DATA, AM_48hrsCNO',2);
[p_anova, tbl, stats] = anova1(DATA)
multcompare(stats)