% Function to Plot LFPBehavior data for all animals in input structure, by
% animal and by period. Subfunction of Analyze_LFPBehavior_CNO_LightDarkCombine_JBmanuscript.m

% INPUT:
%       - LFPBehavior_CombinedData is a structure from Analyze_LFPBehavior_CNO.m (or saved), 
%       which has combined data from all animals listed for the period
%       length listed after CNO injections and baseline period

% PLOTS: 
%       - Behavioral data (percent time, epoch length etc for REM, NREM, AW, QW) after CNO,
%       normalized to baseline
%       - LFP data (delta, theta, alpha power in each behavioral state)
%       after CNO, normalized to baseline period
%       - Plots these things by animal and by CNO period

function stats = Plot_LFPBehavior_CNO_LightDarkCombine_JBmanuscript(STATS)
stats = struct;
stats.PeriodLength = STATS.PeriodLength;
stats.Animals = unique(STATS.ANIMALS_CNOlight);
% plot all the things

%%%%%%%%%% Normalized Plots: CNO Frequency Band Power normalized to Baseline %%%%%%%%%%%%
% BY PERIOD INSTEAD OF BY ANIMAL

REM_color = [0, 0.4470, 0.7410];
NREM_color = [0.4940, 0.1840, 0.5560];
AW_color = [0.8500, 0.3250, 0.0980];
QW_color = [0.6350, 0.0780, 0.1840];
yellow = [0.9290, 0.6940, 0.1250];

normPercents_period = [];

normDurations_period = [];

normTotals_period = [];

normDeltaPeriod = [];
normThetaPeriod = [];
normAlphaPeriod = [];

numPeriods = length(STATS.ANIMALS_CNOlight); % this is assuming that each animal has an equal number of each type of period
PPA = struct;

for pp = 1:numPeriods
    anim = char(STATS.ANIMALS_CNOlight(pp));
    if ~isfield(PPA, anim)
        PPA.(anim) = 0;
    end
    PPA.(anim) = PPA.(anim) + 1;
end

anim_number = 0;
prior = 'none';
prior_num = 0;

for i = 1:numPeriods
    
    anim = char(STATS.ANIMALS_CNOlight(i));
    if ~strcmp(prior, anim)
        anim_number = anim_number+1;
        prior_num = i;
    end
    prior = anim;
    
    normPercents_light = STATS.PercentsPeriod_light_CNO(i,:)./nanmean(STATS.PercentsPeriod_light_baseline(prior_num:prior_num+PPA.(anim)-1,:), 1);
    normPercents_dark = STATS.PercentsPeriod_dark_CNO(i,:)./nanmean(STATS.PercentsPeriod_light_baseline(prior_num:prior_num+PPA.(anim)-1,:), 1);
    normPercents_period(i, :) = nanmean([normPercents_light; normPercents_dark]);
    
    normDurations_light = STATS.DurationsPeriod_light_CNO(i,:)./nanmean(STATS.DurationsPeriod_light_baseline(prior_num:prior_num+PPA.(anim)-1,:),1);
    normDurations_dark = STATS.DurationsPeriod_dark_CNO(i,:)./nanmean(STATS.DurationsPeriod_dark_baseline(prior_num:prior_num+PPA.(anim)-1,:),1);
    normDurations_period(i, :) = nanmean([normDurations_light; normDurations_dark]);
    
    normTotals_light = STATS.totalsPeriod_light_CNO(i,:)./nanmean(STATS.totalsPeriod_light_baseline(prior_num:prior_num+PPA.(anim)-1,:), 1);
    normTotals_dark = STATS.totalsPeriod_dark_CNO(i,:)./nanmean(STATS.totalsPeriod_dark_baseline(prior_num:prior_num+PPA.(anim)-1,:), 1);
    normTotals_period(i, :) = nanmean([normTotals_light; normTotals_dark]);
    
    normDeltaPeriod_light = STATS.DeltaPeriod_light_CNO(i,:)./nanmean(STATS.DeltaPeriod_light_baseline(prior_num:prior_num+PPA.(anim)-1,:), 1);
    normThetaPeriod_light = STATS.ThetaPeriod_light_CNO(i,:)./nanmean(STATS.ThetaPeriod_light_baseline(prior_num:prior_num+PPA.(anim)-1,:), 1);
    normAlphaPeriod_light = STATS.AlphaPeriod_light_CNO(i,:)./nanmean(STATS.AlphaPeriod_light_baseline(prior_num:prior_num+PPA.(anim)-1,:), 1);
    
    normDeltaPeriod_dark = STATS.DeltaPeriod_dark_CNO(i,:)./nanmean(STATS.DeltaPeriod_dark_baseline(prior_num:prior_num+PPA.(anim)-1,:), 1);
    normThetaPeriod_dark = STATS.ThetaPeriod_dark_CNO(i,:)./nanmean(STATS.ThetaPeriod_dark_baseline(prior_num:prior_num+PPA.(anim)-1,:), 1);
    normAlphaPeriod_dark = STATS.AlphaPeriod_dark_CNO(i,:)./nanmean(STATS.AlphaPeriod_dark_baseline(prior_num:prior_num+PPA.(anim)-1,:), 1);
    
    normDeltaPeriod(i, :) = nanmean([normDeltaPeriod_light; normDeltaPeriod_dark]);
    normThetaPeriod(i, :) = nanmean([normThetaPeriod_light; normThetaPeriod_dark]);
    normAlphaPeriod(i, :) = nanmean([normAlphaPeriod_light; normAlphaPeriod_dark]);
end


figure('NumberTitle', 'off', 'Name', 'Percent Time Spent in States by PERIOD')

boxplot (normPercents_period,'BoxStyle', 'filled', 'Whisker', 0, 'Symbol', '', 'Colors', 'k')

ylabel ('CNO/baseline')

%xlabel ('Behavioral State')
set (gca, 'xticklabel', {'REM', 'NREM', 'AW', 'QW'}, 'box', 'off', 'fontsize', 15, 'ylim', [0, 4])
hold on

for p = 1: length(normPercents_period(:,1))
    
plot(1, normPercents_period(p,1), 'o', 'MarkerFaceColor', REM_color, 'MarkerEdgeColor', REM_color, 'MarkerSize', 10)
plot(2, normPercents_period(p,2), 'o', 'MarkerFaceColor', NREM_color, 'MarkerEdgeColor', NREM_color, 'MarkerSize', 10)
plot(3, normPercents_period(p,3), 'o', 'MarkerFaceColor', AW_color, 'MarkerEdgeColor', AW_color, 'MarkerSize', 10)
plot(4, normPercents_period(p,4), 'o', 'MarkerFaceColor', QW_color, 'MarkerEdgeColor', QW_color, 'MarkerSize', 10)

end
plot ([0,4.5], [1,1], 'k:')


figure('NumberTitle', 'off', 'Name', 'Epoch Durations by PERIOD')

boxplot (normDurations_period,'BoxStyle', 'filled', 'Whisker', 0, 'Symbol', '', 'Colors', 'k')

ylabel ('CNO/baseline')
%xlabel ('Behavioral State')
set (gca, 'xticklabel', {'REM', 'NREM', 'AW', 'QW'}, 'box', 'off', 'fontsize', 15, 'ylim', [0, 2])
hold on

for p = 1:length(normDurations_period(:,1))
    
plot(1, normDurations_period(p,1), 'o', 'MarkerFaceColor', REM_color, 'MarkerEdgeColor', REM_color, 'MarkerSize', 10)
plot(2, normDurations_period(p,2), 'o', 'MarkerFaceColor', NREM_color, 'MarkerEdgeColor', NREM_color, 'MarkerSize', 10)
plot(3, normDurations_period(p,3), 'o', 'MarkerFaceColor', AW_color, 'MarkerEdgeColor', AW_color, 'MarkerSize', 10)
plot(4, normDurations_period(p,4), 'o', 'MarkerFaceColor', QW_color, 'MarkerEdgeColor', QW_color, 'MarkerSize', 10)

end
plot ([0,4.5], [1,1], 'k:')

figure('NumberTitle', 'off', 'Name', 'Total Power Changes by PERIOD')

boxplot (normTotals_period,'BoxStyle', 'filled', 'Whisker', 0, 'Symbol', '', 'Colors', 'k')

ylabel ('CNO/baseline')

%xlabel ('Frequency Band')
set (gca, 'xticklabel', {'DELTA (0.5-4 Hz)', 'THETA (6-8 Hz)', 'ALPHA (9-15 Hz)'}, 'box', 'off', 'fontsize', 15, 'ylim', [0.5, 1.5])
hold on

for p = 1:length(normTotals_period(:,1))
    
plot ([1,2,3],[normTotals_period(p,1), normTotals_period(p, 2), normTotals_period(p, 3)], 'o', 'MarkerFaceColor', [.5,.5,.5], 'MarkerEdgeColor', [.5,.5,.5], 'MarkerSize', 10)

end
plot ([0,3.5], [1,1], 'k:')

xlabels = {'REM', 'NREM', 'AW', 'QW'};
figure('NumberTitle', 'off', 'Name', 'LFP-Behavior Analysis: PERIOD')
subplot (1,3,1)
boxplot(normDeltaPeriod, 'BoxStyle', 'filled', 'Whisker', 0, 'Symbol', '', 'Colors', 'k')
ylabel ('CNO/baseline')
title ('DELTA')
%xlabel ('Behavioral State')
set (gca, 'xticklabel', xlabels, 'box', 'off', 'fontsize', 15, 'ylim', [0, 2])
hold on

for i = 1:length(normDeltaPeriod(:,1))

plot(1, normDeltaPeriod(i,1), 'o', 'MarkerFaceColor', REM_color, 'MarkerEdgeColor', REM_color, 'MarkerSize', 10)
plot(2, normDeltaPeriod(i,2), 'o', 'MarkerFaceColor', NREM_color, 'MarkerEdgeColor', NREM_color, 'MarkerSize', 10)
plot(3, normDeltaPeriod(i,3), 'o', 'MarkerFaceColor', AW_color, 'MarkerEdgeColor', AW_color, 'MarkerSize', 10)
plot(4, normDeltaPeriod(i,4), 'o', 'MarkerFaceColor', QW_color, 'MarkerEdgeColor', QW_color, 'MarkerSize', 10)

end
plot ([0,4.5], [1,1], 'k:')

subplot (1,3,2)
boxplot(normThetaPeriod, 'BoxStyle', 'filled', 'Whisker', 0, 'Symbol', '', 'Colors', 'k')
%ylabel ('CNO/baseline')
title ('THETA')
%xlabel ('Behavioral State')
set (gca, 'xticklabel', xlabels, 'box', 'off', 'fontsize', 15, 'ylim', [0, 2])
hold on

for i = 1:length(normThetaPeriod(:,1))

plot(1, normThetaPeriod(i,1), 'o', 'MarkerFaceColor', REM_color, 'MarkerEdgeColor', REM_color, 'MarkerSize', 10)
plot(2, normThetaPeriod(i,2), 'o', 'MarkerFaceColor', NREM_color, 'MarkerEdgeColor', NREM_color, 'MarkerSize', 10)
plot(3, normThetaPeriod(i,3), 'o', 'MarkerFaceColor', AW_color, 'MarkerEdgeColor', AW_color, 'MarkerSize', 10)
plot(4, normThetaPeriod(i,4), 'o', 'MarkerFaceColor', QW_color, 'MarkerEdgeColor', QW_color, 'MarkerSize', 10)

end
plot ([0,4.5], [1,1], 'k:')

subplot (1,3,3)
boxplot(normAlphaPeriod, 'BoxStyle', 'filled', 'Whisker', 0, 'Symbol', '', 'Colors', 'k')
%ylabel ('CNO/baseline')
title ('ALPHA')
%xlabel ('Behavioral State')
set (gca, 'xticklabel', xlabels, 'box', 'off', 'fontsize', 15, 'ylim', [0, 2])
hold on

for i = 1:length(normAlphaPeriod(:,1))

plot(1, normAlphaPeriod(i,1), 'o', 'MarkerFaceColor', REM_color, 'MarkerEdgeColor', REM_color, 'MarkerSize', 10)
plot(2, normAlphaPeriod(i,2), 'o', 'MarkerFaceColor', NREM_color, 'MarkerEdgeColor', NREM_color, 'MarkerSize', 10)
plot(3, normAlphaPeriod(i,3), 'o', 'MarkerFaceColor', AW_color, 'MarkerEdgeColor', AW_color, 'MarkerSize', 10)
plot(4, normAlphaPeriod(i,4), 'o', 'MarkerFaceColor', QW_color, 'MarkerEdgeColor', QW_color, 'MarkerSize', 10)

end
plot ([0,4.5], [1,1], 'k:')

%%% NOTE: for Period Stats order of rows: 2 JB13 (for all light CNO, dark
%%% CNO, light baseline, and dark baseline, in chronological order so first of each baseline is BEFORE CNO injections and second is AFTER), 2
%%% JB14 (same as JB13), JB22 (one light CNO, one dark CNO, two each
%%% baseline, JB26 (two of everything)
%%% Baselines are before the start of the CNO injections for JB22 and JB26)


% STATS for baseline vs CNO periods by Period
for i = 1:4
[stats.DeltaPeriod_dark(1,i), stats.DeltaPeriod_dark(2,i)] = ttest(STATS.DeltaPeriod_dark_baseline(:,i), STATS.DeltaPeriod_dark_CNO(:,i));
end

for i = 1:4
[stats.ThetaPeriod_dark(1,i) ,stats.ThetaPeriod_dark(2,i)] = ttest(STATS.ThetaPeriod_dark_baseline(:,i), STATS.ThetaPeriod_dark_CNO(:,i));
end

for i = 1:4
[stats.AlphaPeriod_dark(1,i), stats.AlphaPeriod_dark(2,i)] = ttest(STATS.AlphaPeriod_dark_baseline(:,i), STATS.AlphaPeriod_dark_CNO(:,i));
end

for i = 1:3
[stats.totalsPeriod_dark(1,i), stats.totalsPeriod_dark(2,i)] = ttest(STATS.totalsPeriod_dark_baseline(:,i), STATS.totalsPeriod_dark_CNO(:,i));
end

for i = 1:3
[stats.totalsPeriod_light(1,i), stats.totalsPeriod_light(2,i)] = ttest(STATS.totalsPeriod_light_baseline(:,i), STATS.totalsPeriod_light_CNO(:,i));
end

for i = 1:4
[stats.DeltaPeriod_light(1,i), stats.DeltaPeriod_light(2,i)] = ttest(STATS.DeltaPeriod_light_baseline(:,i), STATS.DeltaPeriod_light_CNO(:,i));
end

for i = 1:4
[stats.ThetaPeriod_light(1,i), stats.ThetaPeriod_light(2,i)] = ttest(STATS.ThetaPeriod_light_baseline(:,i), STATS.ThetaPeriod_light_CNO(:,i));
end

for i = 1:4
[stats.AlphaPeriod_light(1,i), stats.AlphaPeriod_light(2,i)] = ttest(STATS.AlphaPeriod_light_baseline(:,i), STATS.AlphaPeriod_light_CNO(:,i));
end

for i = 1:4
[stats.PercentsPeriod_light(1,i), stats.PercentsPeriod_light(2,i)] = ttest(STATS.PercentsPeriod_light_baseline(:,i), STATS.PercentsPeriod_light_CNO(:,i));
end

for i = 1:4
[stats.PercentsPeriod_dark(1,i), stats.PercentsPeriod_dark(2,i)] = ttest(STATS.PercentsPeriod_dark_baseline(:,i), STATS.PercentsPeriod_dark_CNO(:,i));
end

for i = 1:4
[stats.DurationsPeriod_light(1,i), stats.DurationsPeriod_light(2,i)] = ttest(STATS.DurationsPeriod_light_baseline(:,i), STATS.DurationsPeriod_light_CNO(:,i));
end

for i = 1:4
[stats.DurationsPeriod_dark(1,i), stats.DurationsPeriod_dark(2,i)] = ttest(STATS.DurationsPeriod_dark_baseline(:,i), STATS.DurationsPeriod_dark_CNO(:,i));
end

%%%%%%%%%% Normalized Plots: CNO Frequency Band Power normalized to Baseline %%%%%%%%%%%%
% BY ANIMAL INSTEAD OF BY PERIOD

REM_color = [0, 0.4470, 0.7410];
NREM_color = [0.4940, 0.1840, 0.5560];
AW_color = [0.8500, 0.3250, 0.0980];
QW_color = [0.6350, 0.0780, 0.1840];
yellow = [0.9290, 0.6940, 0.1250];

Colors = [REM_color; NREM_color; AW_color; QW_color];

normPercents_animal = [];

normDurations_animal = [];

normTotals_animal = [];

normDeltaAnimal = [];
normThetaAnimal = [];
normAlphaAnimal = [];

for i = 1:length(unique(STATS.ANIMALS_CNOlight))
    
    normPercents_light = STATS.PercentsAnimal_light_CNO(i,:)./STATS.PercentsAnimal_light_baseline(i,:);
    normPercents_dark = STATS.PercentsAnimal_dark_CNO(i,:)./STATS.PercentsAnimal_dark_baseline(i,:);
    normPercents_animal(i, :) = nanmean([normPercents_light; normPercents_dark]);
    
    normDurations_light = STATS.DurationsAnimal_light_CNO(i,:)./STATS.DurationsAnimal_light_baseline(i,:);
    normDurations_dark = STATS.DurationsAnimal_dark_CNO(i,:)./STATS.DurationsAnimal_dark_baseline(i,:);
    normDurations_animal(i, :) = nanmean([normDurations_light; normDurations_dark]);
    
    normTotals_light = STATS.totalsAnimal_light_CNO(i,:)./STATS.totalsAnimal_light_baseline(i,:);
    normTotals_dark = STATS.totalsAnimal_dark_CNO(i,:)./STATS.totalsAnimal_dark_baseline(i,:);
    normTotals_animal(i, :) = nanmean([normTotals_light; normTotals_dark]);
    
    normDeltaAnimal_light = STATS.DeltaAnimal_light_CNO(i,:)./nanmean(STATS.DeltaAnimal_light_baseline(i,:));
    normThetaAnimal_light = STATS.ThetaAnimal_light_CNO(i,:)./nanmean(STATS.ThetaAnimal_light_baseline(i,:));
    normAlphaAnimal_light = STATS.AlphaAnimal_light_CNO(i,:)./nanmean(STATS.AlphaAnimal_light_baseline(i,:));
    
    normDeltaAnimal_dark = STATS.DeltaAnimal_dark_CNO(i,:)./nanmean(STATS.DeltaAnimal_dark_baseline(i,:));
    normThetaAnimal_dark = STATS.ThetaAnimal_dark_CNO(i,:)./nanmean(STATS.ThetaAnimal_dark_baseline(i,:));
    normAlphaAnimal_dark = STATS.AlphaAnimal_dark_CNO(i,:)./nanmean(STATS.AlphaAnimal_dark_baseline(i,:));
    
    normDeltaAnimal(i, :) = nanmean([normDeltaAnimal_light; normDeltaAnimal_dark]);
    normThetaAnimal(i, :) = nanmean([normThetaAnimal_light; normThetaAnimal_dark]);
    normAlphaAnimal(i, :) = nanmean([normAlphaAnimal_light; normAlphaAnimal_dark]);
end

figure('NumberTitle', 'off', 'Name', 'Percent Time Spent in State by ANIMAL')

hold on 
plot ([0,4.5], [1,1], 'k:')
p1 = UnivarScatter_ATP(normPercents_animal,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel ('CNO/baseline')
set (gca, 'xticklabel', {'REM', 'NREM', 'AW', 'QW'}, 'box', 'off', 'fontsize', 15, 'ylim', [0, 2])


[p_REM, h] = ttest(normPercents_animal(:,1)-1)
[p_NREM, h] = ttest(normPercents_animal(:,2)-1)
[p_AW, h] = ttest(normPercents_animal(:,3)-1)
[p_QW, h] = ttest(normPercents_animal(:,4)-1)



figure('NumberTitle', 'off', 'Name', 'Absolute Percent Time Spent in State by ANIMAL')
clrs = [1,1,1; 1,0,0];
colormap(clrs)
subplot(1,2,1)
disp ('Light baseline percents:')
nanmean(STATS.PercentsAnimal_light_baseline)
nanstd(STATS.PercentsAnimal_light_baseline)./sqrt(length(STATS.PercentsAnimal_light_baseline(:,1))-1)
disp ('Light CNO percents:')
nanmean(STATS.PercentsAnimal_light_CNO)
nanstd(STATS.PercentsAnimal_light_CNO)./sqrt(length(STATS.PercentsAnimal_light_CNO(:,1))-1)

disp ('Dark baseline percents:')
nanmean(STATS.PercentsAnimal_dark_baseline)
nanstd(STATS.PercentsAnimal_dark_baseline)./sqrt(length(STATS.PercentsAnimal_dark_baseline(:,1))-1)
disp ('Dark CNO percents:')
nanmean(STATS.PercentsAnimal_dark_CNO)
nanstd(STATS.PercentsAnimal_dark_CNO)./sqrt(length(STATS.PercentsAnimal_dark_CNO(:,1))-1)

hold on 
b = bar([nanmean(STATS.PercentsAnimal_light_baseline)', nanmean(STATS.PercentsAnimal_light_CNO)'], 'FaceColor', 'Flat', 'EdgeColor', 'k');

for x = 1:length(STATS.PercentsAnimal_light_baseline(:,1))
    plot([0.85, 1.85, 2.85, 3.85], [STATS.PercentsAnimal_light_baseline(x,1), STATS.PercentsAnimal_light_baseline(x,2), STATS.PercentsAnimal_light_baseline(x,3), STATS.PercentsAnimal_light_baseline(x,4)], 'ko')
    plot([1.15, 2.15, 3.15, 4.15], [STATS.PercentsAnimal_light_CNO(x,1), STATS.PercentsAnimal_light_CNO(x,2), STATS.PercentsAnimal_light_CNO(x,3), STATS.PercentsAnimal_light_CNO(x,4)], 'ko')

end

box off
ylabel('Percent time in state')
title ('LIGHT')
set (gca, 'xticklabel', {'REM', 'NREM', 'AW', 'QW'}, 'box', 'off', 'fontsize', 15)


subplot(1,2,2)

hold on 
b = bar([nanmean(STATS.PercentsAnimal_dark_baseline)', nanmean(STATS.PercentsAnimal_dark_CNO)'], 'FaceColor', 'Flat', 'EdgeColor', 'k');

for x = 1:length(STATS.PercentsAnimal_dark_baseline(:,1))
    plot([0.85, 1.85, 2.85, 3.85], [STATS.PercentsAnimal_dark_baseline(x,1), STATS.PercentsAnimal_dark_baseline(x,2), STATS.PercentsAnimal_dark_baseline(x,3), STATS.PercentsAnimal_dark_baseline(x,4)], 'ko')
    plot([1.15, 2.15, 3.15, 4.15], [STATS.PercentsAnimal_dark_CNO(x,1), STATS.PercentsAnimal_dark_CNO(x,2), STATS.PercentsAnimal_dark_CNO(x,3), STATS.PercentsAnimal_dark_CNO(x,4)], 'ko')

end

box off
ylabel('Percent time in state')
title ('DARK')
set (gca, 'xticklabel', {'REM', 'NREM', 'AW', 'QW'}, 'box', 'off', 'fontsize', 15)










figure('NumberTitle', 'off', 'Name', 'Epoch Durations by ANIMAL')


hold on 
plot ([0,4.5], [1,1], 'k:')
p1 = UnivarScatter_ATP(normDurations_animal,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel ('CNO/baseline')
set (gca, 'xticklabel', {'REM', 'NREM', 'AW', 'QW'}, 'box', 'off', 'fontsize', 15, 'ylim', [0, 2])


figure('NumberTitle', 'off', 'Name', 'Total Power Changes by ANIMAL')

these_colors = [.5,.5,.5; .5,.5,.5; .5,.5,.5];
hold on 
plot ([0,3.5], [1,1], 'k:')
p1 = UnivarScatter_ATP(normTotals_animal,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',these_colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',these_colors);

box off
ylabel ('CNO/baseline')
set (gca, 'xticklabel', {'DELTA (0.5-4 Hz)', 'THETA (6-8 Hz)', 'ALPHA (9-15 Hz)'}, 'box', 'off', 'fontsize', 15, 'ylim', [0.5, 1.5])

[p_delta, h] = ttest(normTotals_animal(:,1) - 1)
[p_theta, h] = ttest(normTotals_animal(:,2) - 1)
[p_alpha, h] = ttest(normTotals_animal(:,3) - 1)

xlabels = {'REM', 'NREM', 'AW', 'QW'};
figure('NumberTitle', 'off', 'Name', 'LFP-Behavior Analysis: BY ANIMAL')
subplot (1,3,1)

hold on 
plot ([0,4.5], [1,1], 'k:')
p1 = UnivarScatter_ATP(normDeltaAnimal,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel ('CNO/baseline')
title ('DELTA')
set (gca, 'xticklabel', {'REM', 'NREM', 'AW', 'QW'}, 'box', 'off', 'fontsize', 15, 'ylim', [0, 2])


subplot (1,3,2)

hold on 
plot ([0,4.5], [1,1], 'k:')
p1 = UnivarScatter_ATP(normThetaAnimal,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
%ylabel ('CNO/baseline')
title ('THETA')
set (gca, 'xticklabel', {'REM', 'NREM', 'AW', 'QW'}, 'box', 'off', 'fontsize', 15, 'ylim', [0, 2])


subplot (1,3,3)

hold on 
plot ([0,4.5], [1,1], 'k:')
p1 = UnivarScatter_ATP(normAlphaAnimal,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
%ylabel ('CNO/baseline')
title ('ALPHA')
set (gca, 'xticklabel', {'REM', 'NREM', 'AW', 'QW'}, 'box', 'off', 'fontsize', 15, 'ylim', [0, 2])



%%% NOTE: for Period Stats order of rows: 2 JB13 (for all light CNO, dark
%%% CNO, light baseline, and dark baseline, in chronological order so first of each baseline is BEFORE CNO injections and second is AFTER), 2
%%% JB14 (same as JB13), JB22 (one light CNO, one dark CNO, two each
%%% baseline, JB26 (two of everything)
%%% Baselines are before the start of the CNO injections for JB22 and JB26)

%%

for i = 1:4
[stats.DeltaAnimal_dark(1,i), stats.DeltaAnimal_dark(2,i)] = ttest(STATS.DeltaAnimal_dark_baseline(:,i), STATS.DeltaAnimal_dark_CNO(:,i));
end

for i = 1:4
[stats.ThetaAnimal_dark(1,i), stats.ThetaAnimal_dark(2,i)] = ttest(STATS.ThetaAnimal_dark_baseline(:,i), STATS.ThetaAnimal_dark_CNO(:,i));
end

for i = 1:4
[stats.AlphaAnimal_dark(1,i), stats.AlphaAnimal_dark(2,i)] = ttest(STATS.AlphaAnimal_dark_baseline(:,i), STATS.AlphaAnimal_dark_CNO(:,i));
end

for i = 1:3
[stats.totalsAnimal_dark(1,i), stats.totalsAnimal_dark(2,i)] = ttest(STATS.totalsAnimal_dark_baseline(:,i), STATS.totalsAnimal_dark_CNO(:,i));
end

for i = 1:3
[stats.totalsAnimal_light(1,i), stats.totalsAnimal_light(2,i)] = ttest(STATS.totalsAnimal_light_baseline(:,i), STATS.totalsAnimal_light_CNO(:,i));
end

for i = 1:4
[stats.DeltaAnimal_light(1,i), stats.DeltaAnimal_light(2,i)] = ttest(STATS.DeltaAnimal_light_baseline(:,i), STATS.DeltaAnimal_light_CNO(:,i));
end

for i = 1:4
[stats.ThetaAnimal_light(1,i), stats.ThetaAnimal_light(2,i)] = ttest(STATS.ThetaAnimal_light_baseline(:,i), STATS.ThetaAnimal_light_CNO(:,i));
end

for i = 1:4
[stats.AlphaAnimal_light(1,i), stats.AlphaAnimal_light(2,i)] = ttest(STATS.AlphaAnimal_light_baseline(:,i), STATS.AlphaAnimal_light_CNO(:,i));
end

for i = 1:4
[stats.PercentsAnimal_light(1,i), stats.PercentsAnimal_light(2,i)] = ttest(STATS.PercentsAnimal_light_baseline(:,i), STATS.PercentsAnimal_light_CNO(:,i));
end

for i = 1:4
[stats.PercentsAnimal_dark(1,i), stats.PercentsAnimal_dark(2,i)] = ttest(STATS.PercentsAnimal_dark_baseline(:,i), STATS.PercentsAnimal_dark_CNO(:,i));
end

for i = 1:4
[stats.DurationsAnimal_light(1,i), stats.DurationsAnimal_light(2,i)] = ttest(STATS.DurationsAnimal_light_baseline(:,i), STATS.DurationsAnimal_light_CNO(:,i));
end

for i = 1:4
[stats.DurationsAnimal_dark(1,i), stats.DuationsAnimal_dark(2,i)] = ttest(STATS.DurationsAnimal_dark_baseline(:,i), STATS.DurationsAnimal_dark_CNO(:,i));
end

end