% Bar graph for % infected ChAT cells in HDB 
% ACh manuscript Figure 6B



JB21 = [12.5, 11.76]; 
JB22 = [66.07, 68.42, 89.71];
JB26 = [94.12, 83.78, 89.36];
JB32 = [78.38, 89.09, 82.5, 88.1, 86.27, 89.66];
JB33 = [84.62, 82.35, 77.5, 90.63];
JB37 = [0];
no_virus = [mean(JB21), JB37];
virus = [mean(JB22), mean(JB26), mean(JB32), mean(JB33)];

figure()
b = bar([mean(no_virus), mean(virus)], 'LineWidth', 3);
hold on

b.FaceColor = 'flat';
b.CData(1,:) = [1,1,1];
b.CData(2,:) = [1,1,1];
ylim([0,100])
set(gca, 'box', 'off', 'FontSize', 15, 'xticklabel', {'CNO only Animals', 'Hm4di + CNO Animals'})
ylabel('Hm4di+ChAT+ / Total ChAT+ * 100')
title ('Hm4Di Viral Infection in HDB/MCPO')
deviation = 0.5;
no_virus_x_jitter = 1 + rand(1,length(no_virus)).*deviation - deviation/2;
virus_x_jitter = 2 + rand(1,length(virus)).*deviation - deviation/2;

for i = 1:length(no_virus)
    plot(no_virus_x_jitter(i), no_virus(i), 'mo', 'MarkerSize', 15, 'MarkerFaceColor', 'm')
end
for i = 1:length(virus)
    plot(virus_x_jitter(i), virus(i), 'bo', 'MarkerSize', 15, 'MarkerFaceColor', 'b')
end
