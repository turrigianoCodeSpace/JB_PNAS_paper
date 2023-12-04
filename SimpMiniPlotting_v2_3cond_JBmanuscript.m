% mini data averaging and plotting
% modified by Juliet Bottorff from Brian Cary and Raul Ramos' script 
% ACh manuscript Figure 2D,E



set(0,'DefaultLineLineWidth',1.5,...
    'DefaultLineMarkerSize',8, ...
    'DefaultAxesLineWidth',2, ...
    'DefaultTextFontName','Arial',...
    'DefaultAxesFontSize',12,...
    'DefaultAxesBox','off',...
    'DefaultAxesFontWeight','Bold');



% From loading ALL_MINI_DATA.mat
data_strct = simple_output;


cond1 = 'Hm4di';
cond2 = 'CNOonly';
cond3 = 'BF_Hm4di';
cond4 = 'BF_CNOonly';

cond1_strct = data_strct.(char(cond1));
cond2_strct = data_strct.(char(cond2));
cond3_strct = data_strct.(char(cond3));
cond4_strct = data_strct.(char(cond4));

% THESE CELLS ALL STARTED DYING/BLOWING UP WITHIN FIRST 20 SWEEPS, SO
% EXCLUDING:
exclude = {'date072221_Cell09',  'date072221_Cell10', 'date072321_Cell04', 'date072321_Cell09', 'date092021_Cell01', 'date092221_Cell09', 'date092221_Cell10', 'date092221_Cell14', 'date092321_Cell04'};

%make this based on formatting of simple output
amp_col = 2;
length_col = 12;

cond1_amp_avgs = [];
cond1_freq_avgs = [];
Hm4di_NUM_minis = {};
count = 0;
Hm4diWAVG = [];
for row = 2:size(cond1_strct,1)
    
    if ismember(cond1_strct{row,1}, exclude)
        continue
    end
    
    % exclude cells based on extreme PP numbers
    %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
        thisCell_PP_data = [];
        for col = 2:12
        
        this_entry = cond1_strct{row,col};
        this_entry_vect = [];
        
        for entry_cell = 1:length(this_entry)
            if ~isempty(this_entry{entry_cell})
            this_entry_vect(entry_cell) = nanmean(this_entry{entry_cell});
            end
        end
        thisCell_PP_data(col-1) = nanmean(this_entry_vect);
        end
    
        if thisCell_PP_data(2) > 20000000 % Ra threshold
            continue
        elseif thisCell_PP_data(3) < 90000000 % Rin threshold
            continue
        elseif thisCell_PP_data(4) < 0.000000000050 % Cp threshold
            continue
        end
        
        % exclude cells that are too noisy
        cell_name = cond1_strct{row,1};
        cell_RMSE = [];
        for t = 1:size(analysis_output.(char(cond1)).(cell_name).mini_data, 1)
        TRACE = analysis_output.(char(cond1)).(cell_name).mini_data{t,2}.DATA;
     
        av = mean(TRACE);

        RMSE = sqrt(mean((TRACE-av).^2));
        cell_RMSE(end+1) = RMSE;
        end
        if mean(cell_RMSE) > 0.0000000000038 
            continue
        end
        
    count = count + 1;
    cell_amps = [];
    for trace = 1:length(cond1_strct{row,amp_col})
        if ~isempty(cond1_strct{row,amp_col}{trace})
        cell_amps = [cell_amps cond1_strct{row,amp_col}{trace}];
        end
    end
    
    total_s = 0;
    cell_num_minis = 0;
    for trace = 1:length(cond1_strct{row,length_col})
        if ~isempty(cond1_strct{row,length_col}{trace})
        tr_length = cond1_strct{row,length_col}{trace};
        num_minis = length(cond1_strct{row,amp_col}{trace});
        
        total_s = total_s + tr_length;
        cell_num_minis = cell_num_minis + num_minis;
        end
    end

    cond1_amp_avgs(end+1) = nanmean(cell_amps);
    cond1_freq_avgs(end+1) = cell_num_minis/total_s;
    
    Hm4di_NUM_minis{count, 1} = cell_name;
    Hm4di_NUM_minis{count, 2} = cell_num_minis;
    Hm4di_NUM_minis{count, 3} = nanmean(cell_amps);
    Hm4di_NUM_minis{count, 4} = mean(cell_RMSE);
    
end

cond2_amp_avgs = [];
cond2_freq_avgs = [];
CNOonly_NUM_minis = {};
count = 0;
CNOonlyWAVG = [];

for row = 2:size(cond2_strct,1)
    
    if ismember(cond2_strct{row,1}, exclude)
        continue
    end

    
        % exclude cells based on extreme PP numbers
        %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
        thisCell_PP_data = [];
        for col = 2:12
        
        this_entry = cond2_strct{row,col};
        this_entry_vect = [];
        
        for entry_cell = 1:length(this_entry)
            if ~isempty(this_entry{entry_cell})
            this_entry_vect(entry_cell) = nanmean(this_entry{entry_cell});
            end
        end
        thisCell_PP_data(col-1) = nanmean(this_entry_vect);
        end
    
        if thisCell_PP_data(2) > 20000000 % Ra
            continue
        elseif thisCell_PP_data(3) < 90000000 % Rin
            continue
        elseif thisCell_PP_data(4) < 0.000000000050 % Cp
            continue
        end
        
        % exclude cells that are too noisy
        cell_name = cond2_strct{row,1};
        cell_RMSE = [];
        for t = 1:size(analysis_output.(char(cond2)).(cell_name).mini_data, 1)
        TRACE = analysis_output.(char(cond2)).(cell_name).mini_data{t,2}.DATA;
    
        av = mean(TRACE);

        RMSE = sqrt(mean((TRACE-av).^2));
        cell_RMSE(end+1) = RMSE;
        end
        if mean(cell_RMSE) > 0.0000000000038
            continue
        end
        
        
        
    count = count + 1;
    cell_amps = [];
    for trace = 1:length(cond2_strct{row,amp_col})
        if ~isempty(cond2_strct{row,amp_col}{trace})
        cell_amps = [cell_amps cond2_strct{row,amp_col}{trace}];
        end
    end
    
    total_s = 0;
    cell_num_minis = 0;
    for trace = 1:length(cond2_strct{row,length_col})
        if ~isempty(cond2_strct{row,length_col}{trace})
        tr_length = cond2_strct{row,length_col}{trace};
        num_minis = length(cond2_strct{row,amp_col}{trace});
        
        total_s = total_s + tr_length;
        cell_num_minis = cell_num_minis + num_minis;
        end
    end
 
    cond2_amp_avgs(end+1) = nanmean(cell_amps);
    cond2_freq_avgs(end+1) = cell_num_minis/total_s;
    
    CNOonly_NUM_minis{count, 1} = cell_name;
    CNOonly_NUM_minis{count, 2} = cell_num_minis;
    CNOonly_NUM_minis{count, 3} = nanmean(cell_amps);
    CNOonly_NUM_minis{count, 4} = mean(cell_RMSE);
    
end

cond3_amp_avgs = [];
cond3_freq_avgs = [];
BF_Hm4di_NUM_minis = {};
count = 0;
BF_Hm4diWAVG = [];
BF_Hm4di_HDBvirus = struct;

for row = 2:size(cond3_strct,1)
    
    if ismember(cond3_strct{row,1}, exclude)
        continue
    end

    
        % exclude cells based on extreme PP numbers
        %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
        thisCell_PP_data = [];
        for col = 2:12
        
        this_entry = cond3_strct{row,col};
        this_entry_vect = [];
        
        for entry_cell = 1:length(this_entry)
            if ~isempty(this_entry{entry_cell})
            this_entry_vect(entry_cell) = nanmean(this_entry{entry_cell});
            end
        end
        thisCell_PP_data(col-1) = nanmean(this_entry_vect);
        end
    
        if thisCell_PP_data(2) > 20000000 
            continue
        elseif thisCell_PP_data(3) < 90000000 
  
            continue
        elseif thisCell_PP_data(4) < 0.000000000050 
            continue
        end
        
        % exclude cells that are too noisy
        cell_name = cond3_strct{row,1};
        cell_RMSE = [];
        for t = 1:size(analysis_output.(char(cond3)).(cell_name).mini_data, 1)
        TRACE = analysis_output.(char(cond3)).(cell_name).mini_data{t,2}.DATA;
      
        av = mean(TRACE);

        RMSE = sqrt(mean((TRACE-av).^2));
        cell_RMSE(end+1) = RMSE;
        end
        if mean(cell_RMSE) > 0.0000000000038
            continue
        end
        
        
        
    count = count + 1;
    cell_amps = [];
    for trace = 1:length(cond3_strct{row,amp_col})
        if ~isempty(cond3_strct{row,amp_col}{trace})
        cell_amps = [cell_amps cond3_strct{row,amp_col}{trace}];
        end
    end
    
    total_s = 0;
    cell_num_minis = 0;
    for trace = 1:length(cond3_strct{row,length_col})
        if ~isempty(cond3_strct{row,length_col}{trace})
        tr_length = cond3_strct{row,length_col}{trace};
        num_minis = length(cond3_strct{row,amp_col}{trace});
        
        total_s = total_s + tr_length;
        cell_num_minis = cell_num_minis + num_minis;
        end
    end

    cond3_amp_avgs(end+1) = nanmean(cell_amps);
    cond3_freq_avgs(end+1) = cell_num_minis/total_s;
    
    BF_Hm4di_NUM_minis{count, 1} = cell_name;
    BF_Hm4di_NUM_minis{count, 2} = cell_num_minis;
    BF_Hm4di_NUM_minis{count, 3} = nanmean(cell_amps);
    BF_Hm4di_NUM_minis{count, 4} = mean(cell_RMSE);
   
end

cond4_amp_avgs = [];
cond4_freq_avgs = [];
BF_CNOonly_NUM_minis = {};
count = 0;
BF_CNOonlyWAVG = [];
BF_CNOonly_HDBvirus = struct;

for row = 2:size(cond4_strct,1)
    
    if ismember(cond4_strct{row,1}, exclude)
        continue
    end

    
        % exclude cells based on extreme PP numbers
        %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
        thisCell_PP_data = [];
        for col = 2:12
        
        this_entry = cond4_strct{row,col};
        this_entry_vect = [];
        
        for entry_cell = 1:length(this_entry)
            if ~isempty(this_entry{entry_cell})
            this_entry_vect(entry_cell) = nanmean(this_entry{entry_cell});
            end
        end
        thisCell_PP_data(col-1) = nanmean(this_entry_vect);
        end
    
        if thisCell_PP_data(2) > 20000000
            continue
        elseif thisCell_PP_data(3) < 90000000
            continue
        elseif thisCell_PP_data(4) < 0.000000000050 
            continue
        end
        
        % exclude cells that are too noisy
        cell_name = cond4_strct{row,1};
        cell_RMSE = [];
        for t = 1:size(analysis_output.(char(cond4)).(cell_name).mini_data, 1)
        TRACE = analysis_output.(char(cond4)).(cell_name).mini_data{t,2}.DATA;

        av = mean(TRACE);

        RMSE = sqrt(mean((TRACE-av).^2));
        cell_RMSE(end+1) = RMSE;
        end
        if mean(cell_RMSE) > 0.0000000000038
            continue
        end
        
        
        
    count = count + 1;
    cell_amps = [];
    for trace = 1:length(cond4_strct{row,amp_col})
        if ~isempty(cond4_strct{row,amp_col}{trace})
        cell_amps = [cell_amps cond4_strct{row,amp_col}{trace}];
        end
    end
    
    total_s = 0;
    cell_num_minis = 0;
    for trace = 1:length(cond4_strct{row,length_col})
        if ~isempty(cond4_strct{row,length_col}{trace})
        tr_length = cond4_strct{row,length_col}{trace};
        num_minis = length(cond4_strct{row,amp_col}{trace});
        
        total_s = total_s + tr_length;
        cell_num_minis = cell_num_minis + num_minis;
        end
    end

    cond4_amp_avgs(end+1) = nanmean(cell_amps);
    cond4_freq_avgs(end+1) = cell_num_minis/total_s;
    
    BF_CNOonly_NUM_minis{count, 1} = cell_name;
    BF_CNOonly_NUM_minis{count, 2} = cell_num_minis;
    BF_CNOonly_NUM_minis{count, 3} = nanmean(cell_amps);
    BF_CNOonly_NUM_minis{count, 4} = mean(cell_RMSE);
    
end


figure_coords = [412 655 252 331];
f1 = figure('Position',figure_coords);

data_toplot = padmat(cond2_amp_avgs',cond1_amp_avgs', 2);
data_toplot = padmat(data_toplot,cond4_amp_avgs', 2);
data_toplot = padmat(data_toplot,cond3_amp_avgs', 2);

Colors = [0 0 0; 1 0 1; 0.5 0.5 0.5; 0.58 0.82 0.98];

p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',2,'MarkerEdgeColor',Colors);

box off
ylabel('mEPSC Amp. (pA)','FontSize',15)
set(gca,'XTickLabel',{'CNO only', 'V1 Hm4di + CNO', 'BF ACh Hm4di + CNO', 'V1 Hm4di + BF ACh Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);

xticks([1 2 3 4])

xlim([0.5 4.5])

ylim([0 15])
title('V1m L23 mini amp avgs','FontSize',15)

p_ranksum_amps_3v4 = ranksum(data_toplot(:,1), data_toplot(:,2))
[p, tbl, stats] = kruskalwallis(data_toplot)
KruskalWallis_amps = multcompare(stats)


figure_coords = [412 655 252 331];

f2 = figure('Position',figure_coords);

data_toplot = padmat(cond2_freq_avgs',cond1_freq_avgs', 2);
data_toplot = padmat(data_toplot, cond4_freq_avgs', 2);
data_toplot = padmat(data_toplot, cond3_freq_avgs', 2);

p1 = UnivarScatter_ATP_JBmanuscript(data_toplot,'BoxType','SEM',...
    'Width',1.5,'Compression',35,'MarkerFaceColor',Colors,...
    'PointSize',40,'StdColor','none','SEMColor',[0 0 0],...
    'Whiskers','lines','WhiskerLineWidth',3,'MarkerEdgeColor',Colors);

box off
ylabel('mEPSC Freq. (Hz)')
set(gca,'XTickLabel',{'CNO only', 'V1 Hm4di + CNO', 'BF ACh Hm4di + CNO', 'V1 Hm4di + BF ACh Hm4di + CNO'},'XTickLabelRotation',45,'FontSize',15);

xticks([1 2 3 4])

xlim([0.5 4.5])

ylim([0 20])

title('V1m L23 mini freq avgs','FontSize',15)

p_ranksum_freqs_1v2 = ranksum(data_toplot(:,1), data_toplot(:,2))
[p, tbl, stats] = kruskalwallis(data_toplot)
kruskalwallis_freqs = multcompare(stats)

