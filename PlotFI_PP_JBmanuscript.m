% PLOT PASSIVE PROPERTIES FROM FI DATA
% ACh manuscript Figure 3,5 Rin plots, plus supplemental info

shRNA = 1; % are you plotting shRNA experiments? 1 = Figure 5; 0 = Figure 3

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

CNOonly_data_mat = NaN(length(FI_CNOonly.filename), 5);
count = 0;
for num = 1:length(FI_CNOonly.filename)
    
    if mean(FI_CNOonly.Ra{1,num}) > 20000000
        continue
    elseif mean(FI_CNOonly.Rin{1,num}) < 90000000 
        continue
    end
    
    count = count + 1;
    %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
    CNOonly_data_mat(num, 2) = nanmean(FI_CNOonly.Ra{1, num});
    CNOonly_data_mat(num, 3) = nanmean(FI_CNOonly.Rin{1, num});
    CNOonly_data_mat(num, 4) = nanmean(FI_CNOonly.Cp{1, num});
    CNOonly_data_mat(num, 5) = nanmean(FI_CNOonly.Vr{1, num});
end


BF_CNOonly_data_mat = NaN(length(FI_BF_CNOonly.filename), 5);
count = 0;
for num = 1:length(FI_BF_CNOonly.filename)
    
    if contains(FI_BF_CNOonly.filename{1,num}, 'date112321') % no good virus expression in BF
        continue
    elseif contains(FI_BF_CNOonly.filename{1,num}, 'date120322_cell_04') % too deep - L4 
        continue
    end
    
    if mean(FI_BF_CNOonly.Ra{1,num}) > 20000000
        continue
    elseif mean(FI_BF_CNOonly.Rin{1,num}) < 90000000 
        continue
    end
    
    count = count + 1;
    %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
    BF_CNOonly_data_mat(num, 2) = nanmean(FI_BF_CNOonly.Ra{1, num});
    BF_CNOonly_data_mat(num, 3) = nanmean(FI_BF_CNOonly.Rin{1, num});
    BF_CNOonly_data_mat(num, 4) = nanmean(FI_BF_CNOonly.Cp{1, num});
    BF_CNOonly_data_mat(num, 5) = nanmean(FI_BF_CNOonly.Vr{1, num});
end


Hm4di_data_mat = NaN(length(FI_Hm4di.filename), 5);
if ~shRNA
count = 0;
for num = 1:length(FI_Hm4di.filename)

    if mean(FI_Hm4di.Ra{1,num}) > 20000000
        continue
    elseif mean(FI_Hm4di.Rin{1,num}) < 90000000 
        continue
    end

    count = count + 1;
    
    Hm4di_data_mat(num, 2) = nanmean(FI_Hm4di.Ra{1, num});
    Hm4di_data_mat(num, 3) = nanmean(FI_Hm4di.Rin{1, num});
    Hm4di_data_mat(num, 4) = nanmean(FI_Hm4di.Cp{1, num});
    Hm4di_data_mat(num, 5) = nanmean(FI_Hm4di.Vr{1, num});
end
end

BF_Hm4di_data_mat = NaN(length(FI_BF_Hm4di.filename), 5);
count = 0;
for num = 1:length(FI_BF_Hm4di.filename)
    
    if contains(FI_BF_Hm4di.filename{1,num}, 'date112321') % no good BF virus expression
        continue
    end

    if mean(FI_BF_Hm4di.Ra{1,num}) > 20000000
        continue
    elseif mean(FI_BF_Hm4di.Rin{1,num}) < 90000000 
        continue
    end

    count = count + 1;
    
    BF_Hm4di_data_mat(num, 2) = nanmean(FI_BF_Hm4di.Ra{1, num});
    BF_Hm4di_data_mat(num, 3) = nanmean(FI_BF_Hm4di.Rin{1, num});
    BF_Hm4di_data_mat(num, 4) = nanmean(FI_BF_Hm4di.Cp{1, num});
    BF_Hm4di_data_mat(num, 5) = nanmean(FI_BF_Hm4di.Vr{1, num});
end


%%%%%
%PassProp plotting
if shRNA
    Colors = [0.5 0.5 0.5; 0.58 0.82 0.98];
else
Colors = [0 0 0; 1 0 1; 0.5 0.5 0.5; 0.58 0.82 0.98];
end

figure_coords = [200,800,250,400];

%%%
%%%
%VC Rin
if ~shRNA
Hm4di_prop_data = Hm4di_data_mat(:,3).*10^-6;
mean_Hm4di_Rin = nanmean(Hm4di_prop_data)
CNOonly_prop_data = CNOonly_data_mat(:,3).*10^-6;
mean_CNOonly_Rin = nanmean(CNOonly_prop_data)
end
BF_CNOonly_prop_data = BF_CNOonly_data_mat(:,3).*10^-6;
mean_BF_CNOonly_Rin = nanmean(BF_CNOonly_prop_data)
BF_Hm4di_prop_data = BF_Hm4di_data_mat(:,3).*10^-6;
mean_BF_Hm4di_Rin = nanmean(BF_Hm4di_prop_data)

if ~shRNA
Hm4di_Rin_sem = nanstd(Hm4di_prop_data)./sqrt(length(Hm4di_prop_data))
CNOonly_Rin_sem = nanstd(CNOonly_prop_data)./sqrt(length(CNOonly_prop_data))
end
BF_Hm4di_Rin_sem = nanstd(BF_Hm4di_prop_data)./sqrt(length(BF_Hm4di_prop_data))
BF_CNOonly_Rin_sem = nanstd(BF_CNOonly_prop_data)./sqrt(length(BF_CNOonly_prop_data))

%fig = figure('Position',figure_coords);
figure()
if ~shRNA
    data_toplot_1 = padmat(CNOonly_prop_data,Hm4di_prop_data, 2);
data_toplot_2 = padmat(BF_CNOonly_prop_data, BF_Hm4di_prop_data, 2);
data_toplot = padmat(data_toplot_1, data_toplot_2, 2);
else
    data_toplot = padmat(BF_CNOonly_prop_data, BF_Hm4di_prop_data, 2);
end
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel('Rin (MegaOhms)')

if shRNA
    set(gca,'XTickLabel',{'m1-3shRNA + CNO', 'm1-3shRNA + Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
else
    set(gca,'XTickLabel',{'Empty Vector + CNO', 'V1 Hm4di + CNO', 'BFACh Hm4di + CNO', 'V1 + BFACh Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
end
if ~shRNA
xticks([1 2 3 4])
else
    xticks([1 2])
end
ylim([0 500])
[p_kruskalwallis_Rin, tbl, stats] = kruskalwallis(data_toplot)
Rin_stats = multcompare(stats)


%%%
%%%
%VC Cp
if ~shRNA
Hm4di_prop_data = Hm4di_data_mat(:,4).*10^12;
mean_Hm4di_Cp = nanmean(Hm4di_prop_data)
CNOonly_prop_data = CNOonly_data_mat(:,4).*10^12;
mean_CNOonly_Cp = nanmean(CNOonly_prop_data)
end
BF_Hm4di_prop_data = BF_Hm4di_data_mat(:,4).*10^12;
mean_BF_Hm4di_Cp = nanmean(BF_Hm4di_prop_data)
BF_CNOonly_prop_data = BF_CNOonly_data_mat(:,4).*10^12;
mean_BF_CNOonly_Cp = nanmean(BF_CNOonly_prop_data)

if ~shRNA
Hm4di_Cp_sem = nanstd(Hm4di_prop_data)./sqrt(length(Hm4di_prop_data))
CNOonly_Cp_sem = nanstd(CNOonly_prop_data)./sqrt(length(CNOonly_prop_data))
end
BF_Hm4di_Cp_sem = nanstd(BF_Hm4di_prop_data)./sqrt(length(BF_Hm4di_prop_data))
BF_CNOonly_Cp_sem = nanstd(BF_CNOonly_prop_data)./sqrt(length(BF_CNOonly_prop_data))


%fig = figure('Position',figure_coords);
figure()
if ~shRNA
    data_toplot_1 = padmat(CNOonly_prop_data,Hm4di_prop_data, 2);
data_toplot_2 = padmat(BF_CNOonly_prop_data, BF_Hm4di_prop_data, 2);
data_toplot = padmat(data_toplot_1, data_toplot_2, 2);
else
    data_toplot = padmat(BF_CNOonly_prop_data, BF_Hm4di_prop_data, 2);
end
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel('Cp (F)')

if shRNA
    set(gca,'XTickLabel',{'m1-3shRNA + CNO', 'm1-3shRNA + Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
else
    set(gca,'XTickLabel',{'Empty Vector + CNO', 'V1 Hm4di + CNO', 'BFACh Hm4di + CNO', 'V1 + BFACh Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
end
if ~shRNA
xticks([1 2 3 4])
else
    xticks([1 2])
end
[p_kruskalwallis_Cp, tbl, stats] = kruskalwallis(data_toplot)
Cp_stats = multcompare(stats)

%%%
%VC Vr
if ~shRNA
Hm4di_prop_data = Hm4di_data_mat(:,5).*10^3;
mean_Hm4di_Vr = nanmean(Hm4di_prop_data)
CNOonly_prop_data = CNOonly_data_mat(:,5).*10^3;
mean_CNOonly_Vr = nanmean(CNOonly_prop_data)
end
BF_Hm4di_prop_data = BF_Hm4di_data_mat(:,5).*10^3;
mean_BF_Hm4di_Vr = nanmean(BF_Hm4di_prop_data)
BF_CNOonly_prop_data = BF_CNOonly_data_mat(:,5).*10^3;
mean_BF_CNOonly_Vr = nanmean(BF_CNOonly_prop_data)

if ~shRNA
Hm4di_Vr_sem = nanstd(Hm4di_prop_data)./sqrt(length(Hm4di_prop_data))
CNOonly_Vr_sem = nanstd(CNOonly_prop_data)./sqrt(length(CNOonly_prop_data))
end
BF_Hm4di_Vr_sem = nanstd(BF_Hm4di_prop_data)./sqrt(length(BF_Hm4di_prop_data))
BF_CNOonly_Vr_sem = nanstd(BF_CNOonly_prop_data)./sqrt(length(BF_CNOonly_prop_data))

%fig = figure('Position',figure_coords);
figure()
if ~shRNA
    data_toplot_1 = padmat(CNOonly_prop_data,Hm4di_prop_data, 2);
data_toplot_2 = padmat(BF_CNOonly_prop_data, BF_Hm4di_prop_data, 2);
data_toplot = padmat(data_toplot_1, data_toplot_2, 2);
else 
    data_toplot = padmat(BF_CNOonly_prop_data, BF_Hm4di_prop_data, 2);
end
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel('V rest (mV)')

if shRNA
    set(gca,'XTickLabel',{'m1-3shRNA + CNO', 'm1-3shRNA + Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
else
    set(gca,'XTickLabel',{'Empty Vector + CNO', 'V1 Hm4di + CNO', 'BFACh Hm4di + CNO', 'V1 + BFACh Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
end
if ~shRNA
xticks([1 2 3 4])
else
    xticks([1 2])
end
[p_kruskalwallis_Vr, tbl, stats] = kruskalwallis(data_toplot)
Vr_stats = multcompare(stats)


%%%
%VC Ra
if ~shRNA
Hm4di_prop_data = Hm4di_data_mat(:,2).*10^-6;
mean_Hm4di_Ra = nanmean(Hm4di_prop_data)
CNOonly_prop_data = CNOonly_data_mat(:,2).*10^-6;
mean_CNOonly_Ra = nanmean(CNOonly_prop_data)
end
BF_Hm4di_prop_data = BF_Hm4di_data_mat(:,2).*10^-6;
mean_BF_Hm4di_Ra = nanmean(BF_Hm4di_prop_data)
BF_CNOonly_prop_data = BF_CNOonly_data_mat(:,2).*10^-6;
mean_BF_CNOonly_Ra = nanmean(BF_CNOonly_prop_data)

if ~shRNA
Hm4di_Ra_sem = nanstd(Hm4di_prop_data)./sqrt(length(Hm4di_prop_data))
CNOonly_Ra_sem = nanstd(CNOonly_prop_data)./sqrt(length(CNOonly_prop_data))
end
BF_Hm4di_Ra_sem = nanstd(BF_Hm4di_prop_data)./sqrt(length(BF_Hm4di_prop_data))
BF_CNOonly_Ra_sem = nanstd(BF_CNOonly_prop_data)./sqrt(length(BF_CNOonly_prop_data))


%fig = figure('Position',figure_coords);
figure()
if ~shRNA
    data_toplot_1 = padmat(CNOonly_prop_data,Hm4di_prop_data, 2);
data_toplot_2 = padmat(BF_CNOonly_prop_data, BF_Hm4di_prop_data, 2);
data_toplot = padmat(data_toplot_1, data_toplot_2, 2);
else
    data_toplot = padmat(BF_CNOonly_prop_data, BF_Hm4di_prop_data, 2);
end
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel('Ra (Ohms)')

if shRNA
    set(gca,'XTickLabel',{'m1-3shRNA + CNO', 'm1-3shRNA + Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
else
    set(gca,'XTickLabel',{'Empty Vector + CNO', 'V1 Hm4di + CNO', 'BFACh Hm4di + CNO', 'V1 + BFACh Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
end
if ~shRNA
xticks([1 2 3 4])
else
    xticks([1 2])
end
[p_kruskalwallis_Ra, tbl, stats] = kruskalwallis(data_toplot)
multcompare(stats)
