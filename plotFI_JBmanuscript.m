% plot aspects of FI curves from .mat files 
% ACh manuscript figures 3 and 5 (minus example traces and Rin - both from other scripts)


shRNA = 0; % are you plotting shRNA experiments? 1 = Figure 5; 0 = Figure 3

animal_ages_names = {'date101921', 'date102021', 'date102121', 'date110121', 'date110321', 'date110421', 'date110521', 'date111821', 'date111921', 'date112321', 'date112421', 'date090822', 'date090922', 'date091322', 'date101422', 'date101522', 'date101622', 'date101822', 'date120222', 'date120322', 'date120622', 'date120722', 'date121422', 'date121522', 'date121922', 'date061723', 'date061823', 'date071023', 'date071123', 'date071223', 'date072423', 'date072523', 'date072623', 'date072723', 'date072823'};
KD_length = [NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 14, 15, 19, 14, 15, 16, 18, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, NaN, 6, 7, 8, 9, 10];

animal_HDBHm4di_names = {'date111821', 'date111921', 'date112421', 'date112921', 'date120121'};
animal_HDBHm4di_quant = [80.14, 125.71, 92.92, 103.98, 93.55];

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

if ~shRNA
Hm4di_IFR_ave = nan(length(FI_Hm4di.filename), length(FI_Hm4di.curr_inj{1,1}));
CNOonly_IFR_ave = nan(length(FI_CNOonly.filename), length(FI_CNOonly.curr_inj{1,1}));
end
BF_Hm4di_IFR_ave = nan(length(FI_BF_Hm4di.filename), length(FI_BF_Hm4di.curr_inj{1,1}));
BF_CNOonly_IFR_ave = nan(length(FI_BF_CNOonly.filename), length(FI_BF_CNOonly.curr_inj{1,1}));

slope_start = 12;
slope_end = 17;
if ~shRNA
Hm4di_IFR_slope = nan(length(FI_Hm4di.filename), 1);
CNOonly_IFR_slope = nan(length(FI_CNOonly.filename), 1);
end
BF_Hm4di_IFR_slope = nan(length(FI_BF_Hm4di.filename), 1);
BF_CNOonly_IFR_slope = nan(length(FI_BF_CNOonly.filename), 1);

if ~shRNA
Hm4di_IFR_AUC = nan(length(FI_Hm4di.filename), 1);
CNOonly_IFR_AUC = nan(length(FI_CNOonly.filename), 1);
end
BF_Hm4di_IFR_AUC = nan(length(FI_BF_Hm4di.filename), 1);
BF_CNOonly_IFR_AUC = nan(length(FI_BF_CNOonly.filename), 1);

if ~shRNA
Hm4di_rheobase = nan(length(FI_Hm4di.filename), 1);
CNOonly_rheobase = nan(length(FI_CNOonly.filename), 1);
end
BF_Hm4di_rheobase = nan(length(FI_BF_Hm4di.filename), 1);
BF_CNOonly_rheobase = nan(length(FI_BF_CNOonly.filename), 1);

if ~shRNA
Hm4di_Rin = nan(length(FI_Hm4di.filename), 1);
CNOonly_Rin = nan(length(FI_CNOonly.filename), 1);
end
BF_Hm4di_Rin = nan(length(FI_BF_Hm4di.filename), 1);
BF_CNOonly_Rin = nan(length(FI_BF_CNOonly.filename), 1);

KDlength = struct;
KDlength.days = [6:20];
KDlength.Hm4di = struct;
KDlength.CNOonly = struct;
KDlength.BF_Hm4di = struct;
KDlength.BF_CNOonly = struct;

BF_Hm4di_HDBvirus = nan(length(FI_BF_Hm4di.filename), 1);
BF_CNOonly_HDBvirus = nan(length(FI_BF_CNOonly.filename), 1);

if ~shRNA
count_Hm4di = 0;
for h = 1:length(FI_Hm4di.filename)


    if mean(FI_Hm4di.Ra{1,h}) > 20000000
        continue
    elseif mean(FI_Hm4di.Rin{1,h}) < 90000000 
        continue

    end

    count_Hm4di = count_Hm4di + 1;
    
    this_name = FI_Hm4di.filename{h}(1:10);
    [Y,L] = ismember(this_name, animal_ages_names);  
    
    Hm4di_IFR_ave(h, :) = FI_Hm4di.IFR_ave{1,h}';
    P = polyfit([slope_start:slope_end],FI_Hm4di.IFR_ave{1,h}(slope_start:slope_end)',1);
    Hm4di_IFR_slope(h) = P(1);
    Hm4di_IFR_AUC(h) = trapz(FI_Hm4di.IFR_ave{1,h}');
    Hm4di_rheobase(h) = FI_Hm4di.rheobase(h);
    Hm4di_Rin(h, :) = mean(FI_Hm4di.Rin{1,h});

    this_KDlength = KD_length(L);
    if ~isnan(this_KDlength)

        if ~isfield(KDlength.Hm4di, this_name)
            KDlength.Hm4di.(this_name) = cell(1,3);
            KDlength.Hm4di.(this_name){1} = this_KDlength;
            KDlength.Hm4di.(this_name){2}(end+1) = Hm4di_IFR_AUC(h);
        else
            KDlength.Hm4di.(this_name){1} = this_KDlength;
            KDlength.Hm4di.(this_name){2}(end+1) = Hm4di_IFR_AUC(h);
        end
    end
    
end

count_CNOonly = 0;
for c = 1:length(FI_CNOonly.filename)

    if contains(FI_CNOonly.filename{1,c}, 'date110121_cell_02') % a few bad sweeps in the middle of current injections (dep up to ~ -40)
        continue
    end
    
    if mean(FI_CNOonly.Ra{1,c}) > 20000000
        continue
    elseif mean(FI_CNOonly.Rin{1,c}) < 90000000
        continue
    end

    count_CNOonly = count_CNOonly + 1;
    
    this_name = FI_CNOonly.filename{c}(1:10);
    [Y,L] = ismember(this_name, animal_ages_names);  
    
    CNOonly_IFR_ave(c, :) = FI_CNOonly.IFR_ave{1,c}';
    P = polyfit([slope_start:slope_end],FI_CNOonly.IFR_ave{1,c}(slope_start:slope_end)',1);
    CNOonly_IFR_slope(c) = P(1);
    CNOonly_IFR_AUC(c) = trapz(FI_CNOonly.IFR_ave{1,c}');
    CNOonly_rheobase(c) = FI_CNOonly.rheobase(c);
    CNOonly_Rin(c, :) = mean(FI_CNOonly.Rin{1,c});

    this_KDlength = KD_length(L);
    if ~isnan(this_KDlength)
   
        if ~isfield(KDlength.CNOonly, this_name)
            KDlength.CNOonly.(this_name) = cell(1,3);
            KDlength.CNOonly.(this_name){1} = this_KDlength;
            KDlength.CNOonly.(this_name){2}(end+1) = CNOonly_IFR_AUC(c);
        else
            KDlength.CNOonly.(this_name){1} = this_KDlength;
            KDlength.CNOonly.(this_name){2}(end+1) = CNOonly_IFR_AUC(c);
        end
    end
    
end
end

count_BF_Hm4di = 0;
for h = 1:length(FI_BF_Hm4di.filename)
    
    if contains(FI_BF_Hm4di.filename{1,h}, 'date112321') % no good virus expression in BF
        continue
    end
    
    if mean(FI_BF_Hm4di.Ra{1,h}) > 20000000
        continue
    elseif mean(FI_BF_Hm4di.Rin{1,h}) < 90000000 
        continue
    end

    count_BF_Hm4di = count_BF_Hm4di + 1;
    
    this_name = FI_BF_Hm4di.filename{h}(1:10);
    [Y,L] = ismember(this_name, animal_ages_names);  
    
    [Yv,Lv] = ismember(this_name, animal_HDBHm4di_names);  
    if Lv == 0
        BF_Hm4di_HDBvirus(h) = NaN;
    else
        BF_Hm4di_HDBvirus(h) = animal_HDBHm4di_quant(Lv);
    end    
    BF_Hm4di_IFR_ave(h, :) = FI_BF_Hm4di.IFR_ave{1,h}';
    P = polyfit([slope_start:slope_end],FI_BF_Hm4di.IFR_ave{1,h}(slope_start:slope_end)',1);
    BF_Hm4di_IFR_slope(h) = P(1);
    BF_Hm4di_IFR_AUC(h) = trapz(FI_BF_Hm4di.IFR_ave{1,h}');
    BF_Hm4di_rheobase(h) = FI_BF_Hm4di.rheobase(h);
    BF_Hm4di_Rin(h, :) = mean(FI_BF_Hm4di.Rin{1,h});

    this_KDlength = KD_length(L);
    if ~isnan(this_KDlength)
    
        if ~isfield(KDlength.BF_Hm4di, this_name)
            KDlength.BF_Hm4di.(this_name) = cell(1,3);
            KDlength.BF_Hm4di.(this_name){1} = this_KDlength;
            KDlength.BF_Hm4di.(this_name){2}(end+1) = BF_Hm4di_IFR_AUC(h);
        else
            KDlength.BF_Hm4di.(this_name){1} = this_KDlength;
            KDlength.BF_Hm4di.(this_name){2}(end+1) = BF_Hm4di_IFR_AUC(h);
        end
    end
    
end

count_BF_CNOonly = 0;
for c = 1:length(FI_BF_CNOonly.filename)
    
    if contains(FI_BF_CNOonly.filename{1,c}, 'date112321') % no good virus expression in BF
        continue
    elseif contains(FI_BF_CNOonly.filename{1,c}, 'date120322_cell_04') % too deep - L4 
        continue
    end
    
    if mean(FI_BF_CNOonly.Ra{1,c}) > 20000000
        continue
    elseif mean(FI_BF_CNOonly.Rin{1,c}) < 90000000
        continue
    end

    count_BF_CNOonly = count_BF_CNOonly + 1;
    
    this_name = FI_BF_CNOonly.filename{c}(1:10);
    [Y,L] = ismember(this_name, animal_ages_names);  
    
    [Yv,Lv] = ismember(this_name, animal_HDBHm4di_names);  
    if Lv == 0
        BF_CNOonly_HDBvirus(c) = NaN;
    else
        BF_CNOonly_HDBvirus(c) = animal_HDBHm4di_quant(Lv);
    end
    
    BF_CNOonly_IFR_ave(c, :) = FI_BF_CNOonly.IFR_ave{1,c}';
    P = polyfit([slope_start:slope_end],FI_BF_CNOonly.IFR_ave{1,c}(slope_start:slope_end)',1);
    BF_CNOonly_IFR_slope(c) = P(1);
    BF_CNOonly_IFR_AUC(c) = trapz(FI_BF_CNOonly.IFR_ave{1,c}');
    BF_CNOonly_rheobase(c) = FI_BF_CNOonly.rheobase(c);
    BF_CNOonly_Rin(c, :) = mean(FI_BF_CNOonly.Rin{1,c});

    this_KDlength = KD_length(L);
    if ~isnan(this_KDlength)

        if ~isfield(KDlength.BF_CNOonly, this_name)
            KDlength.BF_CNOonly.(this_name) = cell(1,3);
            KDlength.BF_CNOonly.(this_name){1} = this_KDlength;
            KDlength.BF_CNOonly.(this_name){2}(end+1) = BF_CNOonly_IFR_AUC(c);
        else
            KDlength.BF_CNOonly.(this_name){1} = this_KDlength;
            KDlength.BF_CNOonly.(this_name){2}(end+1) = BF_CNOonly_IFR_AUC(c);
        end
    end
    
end

if ~shRNA
    Colors = [0 0 0; 1 0 1; 0.5 0.5 0.5; 0.58 0.82 0.98];
else
    Colors = [0.5 0.5 0.5; 0.58 0.82 0.98];
end

i_inj = FI_BF_CNOonly.curr_inj{1,1};

if ~shRNA
mean_Hm4di_IFR = nanmean(Hm4di_IFR_ave, 1);
sem_Hm4di_IFR = nanstd(Hm4di_IFR_ave)./sqrt(count_Hm4di-1);
mean_CNOonly_IFR = nanmean(CNOonly_IFR_ave, 1);
sem_CNOonly_IFR = nanstd(CNOonly_IFR_ave)./sqrt(count_CNOonly-1);
end
mean_BF_Hm4di_IFR = nanmean(BF_Hm4di_IFR_ave, 1);
sem_BF_Hm4di_IFR = nanstd(BF_Hm4di_IFR_ave)./sqrt(count_BF_Hm4di-1);
mean_BF_CNOonly_IFR = nanmean(BF_CNOonly_IFR_ave, 1);
sem_BF_CNOonly_IFR = nanstd(BF_CNOonly_IFR_ave)./sqrt(count_BF_CNOonly-1);

% PLOTS!!

% plot FI curve

figure()
hold on
if ~shRNA
errorbar(i_inj, mean_CNOonly_IFR, sem_CNOonly_IFR, 'Color', Colors(1,:), 'LineWidth', 3)
errorbar(i_inj, mean_Hm4di_IFR, sem_Hm4di_IFR, 'Color', Colors(2,:), 'LineWidth', 3)
errorbar(i_inj, mean_BF_CNOonly_IFR, sem_BF_CNOonly_IFR, 'Color', Colors(3,:), 'LineWidth', 3)
errorbar(i_inj, mean_BF_Hm4di_IFR, sem_BF_Hm4di_IFR, 'Color', Colors(4,:), 'LineWidth', 3)
else
errorbar(i_inj, mean_BF_CNOonly_IFR, sem_BF_CNOonly_IFR, 'Color', Colors(1,:), 'LineWidth', 3)
errorbar(i_inj, mean_BF_Hm4di_IFR, sem_BF_Hm4di_IFR, 'Color', Colors(2,:), 'LineWidth', 3)
end

set(gca, 'box', 'off', 'Fontsize', 20)

if shRNA
    legend('m1-3shRNA + CNO', 'm1-3shRNA + Hm4di + CNO', 'Location', 'NorthWest')
else
    legend('Empty Vector + CNO', 'V1 Hm4di + CNO', 'BFACh Hm4di + CNO', 'V1 + BFACh Hm4di + CNO', 'Location', 'NorthWest')
end
xlabel('Current Injection (pA)')
ylabel('IFR')
    
% plot AUC
figure()
if  ~shRNA
    data_toplot_1 = padmat(CNOonly_IFR_AUC,Hm4di_IFR_AUC, 2);
    data_toplot_2 = padmat(BF_CNOonly_IFR_AUC, BF_Hm4di_IFR_AUC, 2);
    data_toplot = padmat(data_toplot_1, data_toplot_2, 2);
else
    data_toplot = padmat(BF_CNOonly_IFR_AUC, BF_Hm4di_IFR_AUC, 2);
end
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel('IFR AUC')
if shRNA
    set(gca,'XTickLabel',{'m1-3shRNA + CNO', 'm1-3shRNA + Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
else
    set(gca,'XTickLabel',{'Empty vector + CNO', 'V1 Hm4di + CNO', 'BFACh Hm4di + CNO', 'V1 + BFACh Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
end
if ~shRNA
    xticks([1 2 3 4])
else
    xticks([1,2])
end
ylim([0 2000])
[p_kruskalwallis_IFRAUC, tbl, stats] = kruskalwallis(data_toplot)
multcompare(stats)



% plot IFR slope
figure()
if  ~shRNA
    data_toplot_1 = padmat(CNOonly_IFR_slope,Hm4di_IFR_slope, 2);
data_toplot_2 = padmat(BF_CNOonly_IFR_slope, BF_Hm4di_IFR_slope, 2);
data_toplot = padmat(data_toplot_1, data_toplot_2, 2);
else
    data_toplot = padmat(BF_CNOonly_IFR_slope, BF_Hm4di_IFR_slope, 2);
end
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel('IFR Slope')

if shRNA
    set(gca,'XTickLabel',{'m1-3shRNA + CNO', 'm1-3shRNA + Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
else
    set(gca,'XTickLabel',{'Empty Vector + CNO', 'V1 Hm4di + CNO', 'BFACh Hm4di + CNO', 'V1 + BFACh Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
end
if ~shRNA
    xticks([1 2 3 4])
else
    xticks([1,2,])
end
[p_kruskalwallis_IFRslope, tbl, stats] = kruskalwallis(data_toplot)
multcompare(stats)



% plot rheobase

figure()
if ~shRNA
    data_toplot_1 = padmat(CNOonly_rheobase,Hm4di_rheobase, 2);
data_toplot_2 = padmat(BF_CNOonly_rheobase, BF_Hm4di_rheobase,2);
data_toplot = padmat(data_toplot_1, data_toplot_2, 2);
else
    data_toplot = padmat(BF_CNOonly_rheobase, BF_Hm4di_rheobase,2);
end
p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel('Rheobase')

if shRNA
    set(gca,'XTickLabel',{'m1-3shRNA + CNO', 'm1-3shRNA + Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
else
    set(gca,'XTickLabel',{'Empty Vector + CNO', 'V1 Hm4di + CNO', 'BFACh Hm4di + CNO', 'V1 + BFACh Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);
end
ylim([0 300])
if ~shRNA
    xticks([1 2 3 4])
else 
    xticks([1 2])
end
[p_kruskalwallis_rheobase, tbl, stats] = kruskalwallis(data_toplot)
multcompare(stats)


% plot correlation of KD length and IFR AUC DIFFERENCE (M1sh - M1sh/Hm4di) for M1 shRNA expts
if shRNA
figure()
hold on

KD_days_2 = [];
BF_CNOonly_BF_Hm4di = [];


for i = 1:length(fieldnames(KDlength.BF_Hm4di))
    all_names = fieldnames(KDlength.BF_Hm4di);
    this_field = all_names{i};

    if isfield(KDlength.BF_CNOonly, this_field)
        KD_days_2(end+1) = KDlength.BF_Hm4di.(this_field){1};
        BF_CNOonly_BF_Hm4di(end+1) = mean(KDlength.BF_CNOonly.(this_field){2}) - mean(KDlength.BF_Hm4di.(this_field){2});
        plot(KD_days_2(end), BF_CNOonly_BF_Hm4di(end), 'o', 'Color', Colors(2,:), 'MarkerSize', 10, 'LineWidth', 4)
    end
end

if ~isempty(KD_days_2)
[rho_KDlength_2, p_KDlength_2] = corr(KD_days_2', BF_CNOonly_BF_Hm4di');
[poly_KDlength_2] = polyfit(KD_days_2, BF_CNOonly_BF_Hm4di, 1);
plot(KD_days_2, KD_days_2.*poly_KDlength_2(1) + poly_KDlength_2(2), '-', 'Color', Colors(2,:), 'LineWidth', 1)
end

box off
ylabel('(M1sh IFRAUC) - (M1sh/Hm4di IFRAUC)')
xlabel('Days of M1 KD')
set(gca,'FontSize',15, 'box', 'off');
end

% plot correlation of HDB virus quantification and IFR AUC

if ~isempty(BF_Hm4di_HDBvirus(~isnan(BF_Hm4di_HDBvirus)))

figure()
hold on
plot(BF_CNOonly_HDBvirus, BF_CNOonly_IFR_AUC, 'o', 'Color', Colors(3,:), 'MarkerSize', 6, 'LineWidth', 3)
plot(BF_Hm4di_HDBvirus, BF_Hm4di_IFR_AUC, 'o', 'Color', Colors(4,:), 'MarkerSize', 6, 'LineWidth', 3)

[rho_BF_CNOonly_HDBvirus, p_BF_CNOonly_HDBvirus] = corr(BF_CNOonly_HDBvirus(~isnan(BF_CNOonly_HDBvirus)), BF_CNOonly_IFR_AUC(~isnan(BF_CNOonly_IFR_AUC)));
[rho_BF_Hm4di_HDBvirus, p_BF_Hm4di_HDBvirus] = corr(BF_Hm4di_HDBvirus(~isnan(BF_Hm4di_HDBvirus)), BF_Hm4di_IFR_AUC(~isnan(BF_Hm4di_IFR_AUC)));

[poly_BF_CNOonly_HDBvirus] = polyfit(BF_CNOonly_HDBvirus(~isnan(BF_CNOonly_HDBvirus)), BF_CNOonly_IFR_AUC(~isnan(BF_CNOonly_IFR_AUC)), 1);
[poly_BF_Hm4di_HDBvirus] = polyfit(BF_Hm4di_HDBvirus(~isnan(BF_Hm4di_HDBvirus)), BF_Hm4di_IFR_AUC(~isnan(BF_Hm4di_IFR_AUC)), 1);

plot(BF_CNOonly_HDBvirus, BF_CNOonly_HDBvirus.*poly_BF_CNOonly_HDBvirus(1) + poly_BF_CNOonly_HDBvirus(2), '-', 'Color', Colors(3,:), 'LineWidth', 1)
plot(BF_Hm4di_HDBvirus, BF_Hm4di_HDBvirus.*poly_BF_Hm4di_HDBvirus(1) + poly_BF_Hm4di_HDBvirus(2), '-', 'Color', Colors(4,:), 'LineWidth', 1)


box off
ylabel('AUC')
xlabel('HDB Hm4di quantification')
set(gca,'FontSize',15, 'box', 'off');
legend('BFinhib + CNO only', 'BFinhib + Hm4di', 'Location', 'NorthWest')

end
