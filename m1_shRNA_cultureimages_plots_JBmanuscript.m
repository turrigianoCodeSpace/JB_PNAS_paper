% Make cumulative distribution plots of m1 antibody brightness in
% knockdown-transfected/control-transfected cultured neurons
% ACh manuscript figure 4C


CTRL = load('/Users/Juliet/Documents/MATLAB/Turrigiano Lab/CTRL_SH_m1_images.mat');
M1_3 = load('/Users/Juliet/Documents/MATLAB/Turrigiano Lab/m1_3SH_m1_images.mat');

CTRL = CTRL.CTRL_SH_m1_images;
M1_3 = M1_3.m1_3sh_m1_images;

% organize data in structures
CTRL_DATA = struct;
CTRL_DATA.DIV10 = struct;
CTRL_DATA.DIV10.SOMA = struct;

M1_3_DATA = struct;
M1_3_DATA.DIV10 = struct;
M1_3_DATA.DIV10.SOMA = struct;

% Separate out DIV10

% CTRL data
for i = 1:size(CTRL,1)
    
group = CTRL{i,1};
group = group{1}(1:12);

if strcmp(group(11:12), '10') % DIV10
    region = CTRL{i,7};
    if strcmp(region{1}(1), 's') % soma
        if isfield(CTRL_DATA.DIV10.SOMA, ['date_', group(1:6)])
        CTRL_DATA.DIV10.SOMA.(['date_', group(1:6)])(end+1, 1) = CTRL{i,3}{1}; % Area
        CTRL_DATA.DIV10.SOMA.(['date_', group(1:6)])(end, 2) = CTRL{i,4}{1}; % Mean Fluor.
        CTRL_DATA.DIV10.SOMA.(['date_', group(1:6)])(end, 3) = CTRL{i,5}{1}; % Integrated Density
        CTRL_DATA.DIV10.SOMA.(['date_', group(1:6)])(end, 4) = CTRL{i,6}{1}; % Raw Integrated Density
        else
        CTRL_DATA.DIV10.SOMA.(['date_', group(1:6)]) = [];
        CTRL_DATA.DIV10.SOMA.(['date_', group(1:6)])(end+1, 1) = CTRL{i,3}{1}; % Area
        CTRL_DATA.DIV10.SOMA.(['date_', group(1:6)])(end, 2) = CTRL{i,4}{1}; % Mean Fluor.
        CTRL_DATA.DIV10.SOMA.(['date_', group(1:6)])(end, 3) = CTRL{i,5}{1}; % Integrated Density
        CTRL_DATA.DIV10.SOMA.(['date_', group(1:6)])(end, 4) = CTRL{i,6}{1}; % Raw Integrated Density
        end
    end

end

end


% M1_3 data
for i = 1:size(M1_3,1)
    
group = M1_3{i,1};
group = group{1}(1:12);

if strcmp(group(11:12), '10') % DIV10
    region = M1_3{i,7};
    if strcmp(region{1}(1), 's') % soma
        if isfield(M1_3_DATA.DIV10.SOMA, ['date_', group(1:6)])
        M1_3_DATA.DIV10.SOMA.(['date_', group(1:6)])(end+1, 1) = M1_3{i,3}; % Area
        M1_3_DATA.DIV10.SOMA.(['date_', group(1:6)])(end, 2) = M1_3{i,4}; % Mean Fluor.
        M1_3_DATA.DIV10.SOMA.(['date_', group(1:6)])(end, 3) = M1_3{i,5}; % Integrated Density
        M1_3_DATA.DIV10.SOMA.(['date_', group(1:6)])(end, 4) = M1_3{i,6}; % Raw Integrated Density
        else
        M1_3_DATA.DIV10.SOMA.(['date_', group(1:6)]) = [];
        M1_3_DATA.DIV10.SOMA.(['date_', group(1:6)])(end+1, 1) = M1_3{i,3}; % Area
        M1_3_DATA.DIV10.SOMA.(['date_', group(1:6)])(end, 2) = M1_3{i,4}; % Mean Fluor.
        M1_3_DATA.DIV10.SOMA.(['date_', group(1:6)])(end, 3) = M1_3{i,5}; % Integrated Density
        M1_3_DATA.DIV10.SOMA.(['date_', group(1:6)])(end, 4) = M1_3{i,6}; % Raw Integrated Density
        end
    end

end

end

% Normalize to mean CTRL soma or dendrite on that date

%%%%%%%%%%%%%% NORMALIZED VALUES!!! %%%%%%%%%%%%%%%%%%%%%%%%%%%

% norm. mean fluor

CTRL_DIV10_SOMA = [];
CTRL_DIV10_DENDRITE = [];

M1_3_DIV10_SOMA = [];
M1_3_DIV10_DENDRITE = [];

these_fields = fieldnames(CTRL_DATA.DIV10.SOMA);
for d = 1:length(these_fields)
    this_num = length(CTRL_DATA.DIV10.SOMA.(these_fields{d})(:,1));
    THIS_mean = mean(CTRL_DATA.DIV10.SOMA.(these_fields{d})(:,1));
    CTRL_DATA.DIV10.SOMA.(these_fields{d})(1:this_num, 5) = CTRL_DATA.DIV10.SOMA.(these_fields{d})(:,2) ./ THIS_mean; % NORM Mean Fluor
    CTRL_DATA.DIV10.SOMA.(these_fields{d})(1:this_num, 6) = CTRL_DATA.DIV10.SOMA.(these_fields{d})(:,3) ./ THIS_mean; % NORM Integrated Density
    CTRL_DATA.DIV10.SOMA.(these_fields{d})(1:this_num, 7) = CTRL_DATA.DIV10.SOMA.(these_fields{d})(:,4) ./ THIS_mean; % NORM Raw Integrated Density
    CTRL_DIV10_SOMA = [CTRL_DIV10_SOMA, CTRL_DATA.DIV10.SOMA.(these_fields{d})(1:this_num, 5)'];
    
    if isfield(M1_3_DATA.DIV10.SOMA, these_fields{d})
    this_num = length(M1_3_DATA.DIV10.SOMA.(these_fields{d})(:,1));
    M1_3_DATA.DIV10.SOMA.(these_fields{d})(1:this_num, 5) = M1_3_DATA.DIV10.SOMA.(these_fields{d})(:,2) ./ THIS_mean; % NORM Mean Fluor
    M1_3_DATA.DIV10.SOMA.(these_fields{d})(1:this_num, 6) = M1_3_DATA.DIV10.SOMA.(these_fields{d})(:,3) ./ THIS_mean; % NORM Integrated Density
    M1_3_DATA.DIV10.SOMA.(these_fields{d})(1:this_num, 7) = M1_3_DATA.DIV10.SOMA.(these_fields{d})(:,4) ./ THIS_mean; % NORM Raw Integrated Density
    M1_3_DIV10_SOMA = [M1_3_DIV10_SOMA, M1_3_DATA.DIV10.SOMA.(these_fields{d})(1:this_num, 5)'];
    end
end

% Plot cumulative Distributions

% SOMAS
figure()
hold on

[c10, stats_c10] = cdfplot(CTRL_DIV10_SOMA);
c10.Color = 'k';
%c10.LineStyle = '--';
c10.LineWidth = 2;
c10.DisplayName = ['CTRL', ' n = '  num2str(length(CTRL_DIV10_SOMA))];

[m3_10, stats_m3_10] = cdfplot(M1_3_DIV10_SOMA);
m3_10.Color = 'b';
%m3_10.LineStyle = '--';
m3_10.LineWidth = 2;
m3_10.DisplayName = ['M1-3', ' n = '  num2str(length(M1_3_DIV10_SOMA))];
set(gca, 'box', 'off', 'fontsize',  15)
title ('SOMA')
box off
grid off
ylabel ('Cumulative Distribution')
xlabel ('Norm. M1 Fluorescence (A.U.)')
legend('Location', 'Southeast')
[h,p] = kstest2(CTRL_DIV10_SOMA, M1_3_DIV10_SOMA)


