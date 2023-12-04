
% Raw trace plotting
% edited from Brian Cary script by JB
% ACh manuscript Figure 2C

CNOonly_TRACE = ALL_analysis_output.CNOonly.date092221_Cell12.mini_data{8,2}.DATA; % Cell 092221_Cell12 trace 34, t 0.5 to 6.5 seconds
Hm4di_TRACE = ALL_analysis_output.Hm4di.date092321_Cell01.mini_data{5,2}.DATA; % Cell 092321_Cell01 trace 16, t 2.6 to 8.6 seconds
BF_CNOonly_TRACE = ALL_analysis_output.BF_CNOonly.date081821_Cell08.mini_data{2,2}.DATA; % Cell 081821_Cell08 trace 35 t 1.5 to 8.5 seconds
BF_Hm4di_TRACE = ALL_analysis_output.BF_Hm4di.date081821_Cell01.mini_data{4,2}.DATA; % Cell 081821_Cell01 trace 17 t 2.5 to 8.5 seconds

BF_CNOonly = BF_CNOonly_TRACE;
CNOonly = CNOonly_TRACE;
Hm4di = Hm4di_TRACE;
BF_Hm4di = BF_Hm4di_TRACE;


BF_CNOonly = BF_CNOonly - nanmean(BF_CNOonly);
CNOonly = CNOonly - nanmean(CNOonly);
Hm4di = Hm4di - nanmean(Hm4di);

samp_rate = 10e3;
time = 1/samp_rate:1/samp_rate:length(CNOonly)/samp_rate;

colors = [0 0 0; 1 0 1; 0.5 0.5 0.5; 0.58 0.82 0.98];
figure_coords = [669 770 728 208];


f1 = figure('Position',figure_coords);
plot(time,sgolayfilt(BF_CNOonly,2,7),'Color',colors(3,:));
xlim_left = 1.5;
xlim_right = 7.5;
xlim([xlim_left xlim_right])
ylim([-31e-12, 10e-12])
box off;axis off;
hold on
plot([xlim_right-1.6 xlim_right-0.6], [-30e-12 -30e-12], 'k', 'LineWidth', 4)
plot([xlim_right-0.6 xlim_right-0.6], [-20e-12 -30.5e-12], 'k', 'LineWidth', 4)


f2 = figure('Position',figure_coords);
plot(time,sgolayfilt(CNOonly,2,7),'Color',colors(1,:));
xlim_left = 2;
xlim_right = 8;
xlim([xlim_left xlim_right])
ylim([-31e-12, 10e-12])
box off;axis off;
hold on
plot([xlim_right-1.6 xlim_right-0.6], [-30e-12 -30e-12], 'k', 'LineWidth', 4)
plot([xlim_right-0.6 xlim_right-0.6], [-20e-12 -30.5e-12], 'k', 'LineWidth', 4)

f3 = figure('Position',figure_coords);
plot(time,sgolayfilt(Hm4di,2,7),'Color',colors(2,:));
xlim_left = 2.5;
xlim_right = 8.5;
xlim([xlim_left xlim_right])
ylim([-31e-12, 10e-12])
box off;axis off; hold on

plot([xlim_right-1.6 xlim_right-0.6], [-30e-12 -30e-12], 'k', 'LineWidth', 4)
plot([xlim_right-0.6 xlim_right-0.6], [-20e-12 -30.5e-12], 'k', 'LineWidth', 4)

f4 = figure('Position',figure_coords);
plot(time,sgolayfilt(BF_Hm4di,2,7),'Color',colors(4,:));
xlim_left = 2.5;
xlim_right = 8.5;
xlim([xlim_left xlim_right])
ylim([-31e-12, 10e-12])
box off;axis off; hold on

plot([xlim_right-1.6 xlim_right-0.6], [-30e-12 -30e-12], 'k', 'LineWidth', 4)
plot([xlim_right-0.6 xlim_right-0.6], [-20e-12 -30.5e-12], 'k', 'LineWidth', 4)




