%
% Extracting data from data structures and plotting Passive Properties of
% minis
% JB edits from BC
% ACh manuscript supplemental info
%%%%%%%%%%


% THESE CELLS ALL STARTED DYING/BLOWING UP WITHIN FIRST 20 SWEEPS, SO
% EXCLUDING:
exclude = {'date072221_Cell09',  'date072221_Cell10', 'date072321_Cell04', 'date072321_Cell09', 'date092021_Cell01', 'date092221_Cell09', 'date092221_Cell10', 'date092221_Cell14', 'date092321_Cell04'};

% From loading ALL_MINI_DATA.mat
data_strct = simple_output;

CNOonly_strct_fc = data_strct.CNOonly;
Hm4di_strct_fc = data_strct.Hm4di;
BF_CNOonly_strct = data_strct.BF_CNOonly;
BF_Hm4di_strct = data_strct.BF_Hm4di;

CNOonly_data_mat = [];
count = 0;
for row = 2:size(CNOonly_strct_fc,1)
    
    if ismember(CNOonly_strct_fc{row,1}, exclude)
        continue
    end
        % exclude cells based on extreme PP numbers
        %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
        thisCell_PP_data = [];
        for col = 2:12
        
        this_entry = CNOonly_strct_fc{row,col};
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
        cell_name = CNOonly_strct_fc{row,1};
        cell_RMSE = [];
        for t = 1:size(analysis_output.CNOonly.(cell_name).mini_data, 1)
        TRACE = analysis_output.CNOonly.(cell_name).mini_data{t,2}.DATA;
        av = mean(TRACE);

        RMSE = sqrt(mean((TRACE-av).^2));
        cell_RMSE(end+1) = RMSE;
        end
        
        if mean(cell_RMSE) > 0.0000000000038
            continue
        end
        
    
    count = count + 1;
    for col = 2:12
        
        this_entry = CNOonly_strct_fc{row,col};
        this_entry_vect = [];
        
        for entry_cell = 1:length(this_entry)
            if ~isempty(this_entry{entry_cell})
            this_entry_vect(entry_cell) = nanmean(this_entry{entry_cell});
            end
        end
        
        CNOonly_data_mat(count,col-1) = nanmean(this_entry_vect);
        
    end
end


BF_CNOonly_data_mat = [];
count = 0;
for row = 2:size(BF_CNOonly_strct,1)
    
    if ismember(BF_CNOonly_strct{row,1}, exclude)
        continue
    end
    
        % exclude cells based on extreme PP numbers
        %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
        thisCell_PP_data = [];
        for col = 2:12
        
        this_entry = BF_CNOonly_strct{row,col};
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
        cell_name = BF_CNOonly_strct{row,1};
        cell_RMSE = [];
        for t = 1:size(analysis_output.BF_CNOonly.(cell_name).mini_data, 1)
        TRACE = analysis_output.BF_CNOonly.(cell_name).mini_data{t,2}.DATA;
        av = mean(TRACE);

        RMSE = sqrt(mean((TRACE-av).^2));
        cell_RMSE(end+1) = RMSE;
        end
        
        if mean(cell_RMSE) > 0.0000000000038
            continue
        end
        
    
    count = count + 1;
    for col = 2:12
        
        this_entry = BF_CNOonly_strct{row,col};
        this_entry_vect = [];
        
        for entry_cell = 1:length(this_entry)
            if ~isempty(this_entry{entry_cell})
            this_entry_vect(entry_cell) = nanmean(this_entry{entry_cell});
            end
        end
        
        BF_CNOonly_data_mat(count,col-1) = nanmean(this_entry_vect);
        
    end
end

Hm4di_data_mat = [];
count = 0;
for row = 2:size(Hm4di_strct_fc,1)

    if ismember(Hm4di_strct_fc{row,1}, exclude)
        continue
    end
    
        % exclude cells based on extreme PP numbers
        %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
        thisCell_PP_data = [];
        for col = 2:12
        
        this_entry = Hm4di_strct_fc{row,col};
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
        cell_name = Hm4di_strct_fc{row,1};
        cell_RMSE = [];
        for t = 1:size(analysis_output.Hm4di.(cell_name).mini_data, 1)
        TRACE = analysis_output.Hm4di.(cell_name).mini_data{t,2}.DATA;
        av = mean(TRACE);

        RMSE = sqrt(mean((TRACE-av).^2));
        cell_RMSE(end+1) = RMSE;
        end
        
        if mean(cell_RMSE) > 0.0000000000038
            continue
        end
        
        
    count = count + 1;
    for col = 2:12
        
        this_entry = Hm4di_strct_fc{row,col};
        this_entry_vect = [];
        
        for entry_cell = 1:length(this_entry)
            if ~isempty(this_entry{entry_cell})
            this_entry_vect(entry_cell) = nanmean(this_entry{entry_cell});
            end
        end
        
        Hm4di_data_mat(count,col-1) = nanmean(this_entry_vect);

    end
end

BF_Hm4di_data_mat = [];
count = 0;
for row = 2:size(BF_Hm4di_strct,1)

    if ismember(BF_Hm4di_strct{row,1}, exclude)
        continue
    end
    
        % exclude cells based on extreme PP numbers
        %[VC_amps, VC_Ra, VC_Rin, VC_Cp, VC_Vr]
        thisCell_PP_data = [];
        for col = 2:12
        
        this_entry = BF_Hm4di_strct{row,col};
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
        cell_name = BF_Hm4di_strct{row,1};
        cell_RMSE = [];
        for t = 1:size(analysis_output.BF_Hm4di.(cell_name).mini_data, 1)
        TRACE = analysis_output.BF_Hm4di.(cell_name).mini_data{t,2}.DATA;
        av = mean(TRACE);

        RMSE = sqrt(mean((TRACE-av).^2));
        cell_RMSE(end+1) = RMSE;
        end
        
        if mean(cell_RMSE) > 0.0000000000038
            continue
        end
        
        
    count = count + 1;
    for col = 2:12
        
        this_entry = BF_Hm4di_strct{row,col};
        this_entry_vect = [];
        
        for entry_cell = 1:length(this_entry)
            if ~isempty(this_entry{entry_cell})
            this_entry_vect(entry_cell) = nanmean(this_entry{entry_cell});
            end
        end
        
        BF_Hm4di_data_mat(count,col-1) = nanmean(this_entry_vect);

    end
end

%%%%%
%PassProp plotting
Colors = [0 0 0; 1 0 1; 0.5 0.5 0.5; 0.58 0.82 0.98];

figure_coords = [200,800,250,400];

%%%
%%%
%VC Rin
Hm4di_prop_data = Hm4di_data_mat(:,3).*10^-6;
mean_Hm4di_Rin = nanmean(Hm4di_prop_data)
CNOonly_prop_data = CNOonly_data_mat(:,3).*10^-6;
mean_CNOonly_Rin = nanmean(CNOonly_prop_data)
BF_CNOonly_prop_data = BF_CNOonly_data_mat(:,3).*10^-6;
mean_BF_CNOonly_Rin = nanmean(BF_CNOonly_prop_data)
BF_Hm4di_prop_data = BF_Hm4di_data_mat(:,3).*10^-6;
mean_BF_Hm4di_Rin = nanmean(BF_Hm4di_prop_data)


Hm4di_prop_data_sem = nanstd(Hm4di_prop_data)./sqrt(length(Hm4di_prop_data));
Hm4di_SEM_Rin = Hm4di_prop_data_sem
CNOonly_prop_data_sem = nanstd(CNOonly_prop_data)./sqrt(length(CNOonly_prop_data));
CNOonly_SEM_Rin = CNOonly_prop_data_sem
BF_Hm4di_prop_data_sem = nanstd(BF_Hm4di_prop_data)./sqrt(length(BF_Hm4di_prop_data));
BF_Hm4di_SEM_Rin = BF_Hm4di_prop_data_sem
BF_CNOonly_prop_data_sem = nanstd(BF_CNOonly_prop_data)./sqrt(length(BF_CNOonly_prop_data));
BF_CNOonly_SEM_Rin = BF_CNOonly_prop_data_sem

deviation = 0.35;
Hm4di_x_jitter = 2 + rand(1,length(Hm4di_prop_data)).*deviation - deviation/2;
CNOonly_x_jitter = 1 + rand(1,length(CNOonly_prop_data)).*deviation - deviation/2;
BF_Hm4di_x_jitter = 4 + rand(1,length(BF_Hm4di_prop_data)).*deviation - deviation/2;
BF_CNOonly_x_jitter = 3 + rand(1,length(BF_CNOonly_prop_data)).*deviation - deviation/2;

fig = figure('Position',figure_coords);
plot(Hm4di_x_jitter, Hm4di_prop_data, 'o', 'Color', Colors(2,:));
hold on;
plot(CNOonly_x_jitter, CNOonly_prop_data, 'o', 'Color', Colors(1,:));
plot(BF_CNOonly_x_jitter, BF_CNOonly_prop_data, 'o', 'Color', Colors(3,:));
plot(BF_Hm4di_x_jitter, BF_Hm4di_prop_data, 'o', 'Color', Colors(4,:));

bar(2,nanmean(Hm4di_prop_data),'EdgeColor', Colors(2,:),'LineWidth',2,'FaceAlpha',0.0)
bar(1,nanmean(CNOonly_prop_data),'EdgeColor', Colors(1,:),'LineWidth',2,'FaceAlpha',0.0)
bar(4,nanmean(BF_Hm4di_prop_data),'EdgeColor', Colors(4,:),'LineWidth',2,'FaceAlpha',0.0)
bar(3,nanmean(BF_CNOonly_prop_data),'EdgeColor', Colors(3,:),'LineWidth',2,'FaceAlpha',0.0)

errorbar(2, nanmean(Hm4di_prop_data), Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(1, nanmean(CNOonly_prop_data), CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(4, nanmean(BF_Hm4di_prop_data), BF_Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(3, nanmean(BF_CNOonly_prop_data), BF_CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
box off
ylabel('VC Rin (Ohms)')
set(gca,'XTickLabel',{'CNO only', 'Hm4di + CNO', 'BFinhib + CNO only', 'BFinhib + Hm4di'}, 'FontSize', 15)
xticks([1 2 3 4])
xtickangle(45)

data = padmat(CNOonly_prop_data,Hm4di_prop_data, 2);
data = padmat(data, BF_CNOonly_prop_data, 2);
data = padmat(data, BF_Hm4di_prop_data, 2);
[p, tbl, stats] = kruskalwallis(data)
kruskalwallis_Rin = multcompare(stats)

%%%
%%%
%VC Cp
Hm4di_prop_data = Hm4di_data_mat(:,4).*10^12;
mean_Hm4di_Cp = nanmean(Hm4di_prop_data)
CNOonly_prop_data = CNOonly_data_mat(:,4).*10^12;
mean_CNOonly_Cp = nanmean(CNOonly_prop_data)
BF_Hm4di_prop_data = BF_Hm4di_data_mat(:,4).*10^12;
mean_BF_Hm4di_Cp = nanmean(BF_Hm4di_prop_data)
BF_CNOonly_prop_data = BF_CNOonly_data_mat(:,4).*10^12;
mean_BF_CNOonly_Cp = nanmean(BF_CNOonly_prop_data)

Hm4di_prop_data_sem = nanstd(Hm4di_prop_data)./sqrt(length(Hm4di_prop_data));
Hm4di_SEM_Cp = Hm4di_prop_data_sem
CNOonly_prop_data_sem = nanstd(CNOonly_prop_data)./sqrt(length(CNOonly_prop_data));
CNOonly_SEM_Cp = CNOonly_prop_data_sem
BF_Hm4di_prop_data_sem = nanstd(BF_Hm4di_prop_data)./sqrt(length(BF_Hm4di_prop_data));
BF_Hm4di_SEM_Cp = BF_Hm4di_prop_data_sem
BF_CNOonly_prop_data_sem = nanstd(BF_CNOonly_prop_data)./sqrt(length(BF_CNOonly_prop_data));
BF_CNOonly_SEM_Cp = BF_CNOonly_prop_data_sem


deviation = 0.35;
Hm4di_x_jitter = 2 + rand(1,length(Hm4di_prop_data)).*deviation - deviation/2;
CNOonly_x_jitter = 1 + rand(1,length(CNOonly_prop_data)).*deviation - deviation/2;
BF_Hm4di_x_jitter = 4 + rand(1,length(BF_Hm4di_prop_data)).*deviation - deviation/2;
BF_CNOonly_x_jitter = 3 + rand(1,length(BF_CNOonly_prop_data)).*deviation - deviation/2;

fig = figure('Position',figure_coords);
plot(Hm4di_x_jitter, Hm4di_prop_data, 'o', 'Color', Colors(2,:));
hold on;
plot(CNOonly_x_jitter, CNOonly_prop_data, 'o', 'Color', Colors(1,:));
plot(BF_CNOonly_x_jitter, BF_CNOonly_prop_data, 'o', 'Color', Colors(3,:));
plot(BF_Hm4di_x_jitter, BF_Hm4di_prop_data, 'o', 'Color', Colors(4,:));

bar(2,nanmean(Hm4di_prop_data),'EdgeColor',Colors(2,:),'LineWidth',2,'FaceAlpha',0.0)
bar(1,nanmean(CNOonly_prop_data),'EdgeColor',Colors(1,:),'LineWidth',2,'FaceAlpha',0.0)
bar(4,nanmean(BF_Hm4di_prop_data),'EdgeColor',Colors(4,:),'LineWidth',2,'FaceAlpha',0.0)
bar(3,nanmean(BF_CNOonly_prop_data),'EdgeColor',Colors(3,:),'LineWidth',2,'FaceAlpha',0.0)

errorbar(2, nanmean(Hm4di_prop_data), Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(1, nanmean(CNOonly_prop_data), CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(4, nanmean(BF_Hm4di_prop_data), BF_Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(3, nanmean(BF_CNOonly_prop_data), BF_CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
box off
ylabel('VC Cp (F)')
set(gca,'XTickLabel',{'CNO only', 'Hm4di + CNO', 'BFinhib + CNO only', 'BFinhib + Hm4di'}, 'FontSize', 15)
xticks([1 2 3 4])
xtickangle(45)


data = padmat(CNOonly_prop_data,Hm4di_prop_data, 2);
data = padmat(data, BF_CNOonly_prop_data, 2);
data = padmat(data, BF_Hm4di_prop_data, 2);
[p, tbl, stats] = kruskalwallis(data)
kruskalwallis_Cp = multcompare(stats)

%%%
%VC Vr
Hm4di_prop_data = Hm4di_data_mat(:,5).*10^3;
mean_Hm4di_Vr = nanmean(Hm4di_prop_data)
CNOonly_prop_data = CNOonly_data_mat(:,5).*10^3;
mean_CNOonly_Vr = nanmean(CNOonly_prop_data)
BF_Hm4di_prop_data = BF_Hm4di_data_mat(:,5).*10^3;
mean_BF_Hm4di_Vr = nanmean(BF_Hm4di_prop_data)
BF_CNOonly_prop_data = BF_CNOonly_data_mat(:,5).*10^3;
mean_BF_CNOonly_Vr = nanmean(BF_CNOonly_prop_data)


Hm4di_prop_data_sem = nanstd(Hm4di_prop_data)./sqrt(length(Hm4di_prop_data));
Hm4di_SEM_Vr = Hm4di_prop_data_sem
CNOonly_prop_data_sem = nanstd(CNOonly_prop_data)./sqrt(length(CNOonly_prop_data));
CNOonly_SEM_Vr = CNOonly_prop_data_sem
BF_Hm4di_prop_data_sem = nanstd(BF_Hm4di_prop_data)./sqrt(length(BF_Hm4di_prop_data));
BF_Hm4di_SEM_Vr = BF_Hm4di_prop_data_sem
BF_CNOonly_prop_data_sem = nanstd(BF_CNOonly_prop_data)./sqrt(length(BF_CNOonly_prop_data));
BF_CNOonly_SEM_Vr = BF_CNOonly_prop_data_sem

deviation = 0.35;
Hm4di_x_jitter = 2 + rand(1,length(Hm4di_prop_data)).*deviation - deviation/2;
CNOonly_x_jitter = 1 + rand(1,length(CNOonly_prop_data)).*deviation - deviation/2;
BF_Hm4di_x_jitter = 4 + rand(1,length(BF_Hm4di_prop_data)).*deviation - deviation/2;
BF_CNOonly_x_jitter = 3 + rand(1,length(BF_CNOonly_prop_data)).*deviation - deviation/2;

fig = figure('Position',figure_coords);
plot(Hm4di_x_jitter, Hm4di_prop_data, 'o', 'Color', Colors(2,:));
hold on;
plot(CNOonly_x_jitter, CNOonly_prop_data, 'o', 'Color', Colors(1,:));
plot(BF_CNOonly_x_jitter, BF_CNOonly_prop_data, 'o', 'Color', Colors(3,:));
plot(BF_Hm4di_x_jitter, BF_Hm4di_prop_data, 'o', 'Color', Colors(4,:));

bar(2,nanmean(Hm4di_prop_data),'EdgeColor', Colors(2,:),'LineWidth',2,'FaceAlpha',0.0)
bar(1,nanmean(CNOonly_prop_data),'EdgeColor', Colors(1,:),'LineWidth',2,'FaceAlpha',0.0)
bar(4,nanmean(BF_Hm4di_prop_data),'EdgeColor', Colors(4,:),'LineWidth',2,'FaceAlpha',0.0)
bar(3,nanmean(BF_CNOonly_prop_data),'EdgeColor', Colors(3,:),'LineWidth',2,'FaceAlpha',0.0)

errorbar(2, nanmean(Hm4di_prop_data), Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(1, nanmean(CNOonly_prop_data), CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(4, nanmean(BF_Hm4di_prop_data), BF_Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(3, nanmean(BF_CNOonly_prop_data), BF_CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
box off
ylabel('VC Vr (mV)')
set(gca,'XTickLabel',{'CNO only','Hm4di + CNO', 'BFinhib + CNO only', 'BFinhib + Hm4di'}, 'FontSize', 15)
xticks([1 2 3 4])
xtickangle(45)


data = padmat(CNOonly_prop_data,Hm4di_prop_data, 2);
data = padmat(data, BF_CNOonly_prop_data, 2);
data = padmat(data, BF_Hm4di_prop_data, 2);
[p, tbl, stats] = kruskalwallis(data)
kruskalwallis_Vr = multcompare(stats)

%%%
%VC Ra
Hm4di_prop_data = Hm4di_data_mat(:,2).*10^-6;
mean_Hm4di_Ra = nanmean(Hm4di_prop_data)
CNOonly_prop_data = CNOonly_data_mat(:,2).*10^-6;
mean_CNOonly_Ra = nanmean(CNOonly_prop_data)
BF_Hm4di_prop_data = BF_Hm4di_data_mat(:,2).*10^-6;
mean_BF_Hm4di_Ra = nanmean(BF_Hm4di_prop_data)
BF_CNOonly_prop_data = BF_CNOonly_data_mat(:,2).*10^-6;
mean_BF_CNOonly_Ra = nanmean(BF_CNOonly_prop_data)

Hm4di_prop_data_sem = nanstd(Hm4di_prop_data)./sqrt(length(Hm4di_prop_data));
Hm4di_SEM_Ra = Hm4di_prop_data_sem
CNOonly_prop_data_sem = nanstd(CNOonly_prop_data)./sqrt(length(CNOonly_prop_data));
CNOonly_SEM_Ra = CNOonly_prop_data_sem
BF_Hm4di_prop_data_sem = nanstd(BF_Hm4di_prop_data)./sqrt(length(BF_Hm4di_prop_data));
BF_Hm4di_SEM_Ra = BF_Hm4di_prop_data_sem
BF_CNOonly_prop_data_sem = nanstd(BF_CNOonly_prop_data)./sqrt(length(BF_CNOonly_prop_data));
BF_CNOonly_SEM_Ra = BF_CNOonly_prop_data_sem

deviation = 0.35;
Hm4di_x_jitter = 2 + rand(1,length(Hm4di_prop_data)).*deviation - deviation/2;
CNOonly_x_jitter = 1 + rand(1,length(CNOonly_prop_data)).*deviation - deviation/2;
BF_Hm4di_x_jitter = 4 + rand(1,length(BF_Hm4di_prop_data)).*deviation - deviation/2;
BF_CNOonly_x_jitter = 3 + rand(1,length(BF_CNOonly_prop_data)).*deviation - deviation/2;

fig = figure('Position',figure_coords);
plot(Hm4di_x_jitter, Hm4di_prop_data, 'o', 'Color', Colors(2,:));
hold on;
plot(CNOonly_x_jitter, CNOonly_prop_data, 'o', 'Color', Colors(1,:));
plot(BF_CNOonly_x_jitter, BF_CNOonly_prop_data, 'o', 'Color', Colors(3,:));
plot(BF_Hm4di_x_jitter, BF_Hm4di_prop_data, 'o', 'Color', Colors(4,:));

bar(2,nanmean(Hm4di_prop_data),'EdgeColor', Colors(2,:),'LineWidth',2,'FaceAlpha',0.0)
bar(1,nanmean(CNOonly_prop_data),'EdgeColor', Colors(1,:),'LineWidth',2,'FaceAlpha',0.0)
bar(4,nanmean(BF_Hm4di_prop_data),'EdgeColor', Colors(4,:),'LineWidth',2,'FaceAlpha',0.0)
bar(3,nanmean(BF_CNOonly_prop_data),'EdgeColor', Colors(3,:),'LineWidth',2,'FaceAlpha',0.0)

errorbar(2, nanmean(Hm4di_prop_data), Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(1, nanmean(CNOonly_prop_data), CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(4, nanmean(BF_Hm4di_prop_data), BF_Hm4di_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
errorbar(3, nanmean(BF_CNOonly_prop_data), BF_CNOonly_prop_data_sem, 'k','MarkerSize', 15, 'LineWidth', 3.5,'CapSize',20);
box off
ylabel('VC Ra (Ohms)')
set(gca,'XTickLabel',{'CNO only', 'Hm4di + CNO', 'BFinhib + CNO only', 'BFinhib + Hm4di'}, 'FontSize', 15)
xticks([1 2 3 4])
xtickangle(45)

data = padmat(CNOonly_prop_data,Hm4di_prop_data, 2);
data = padmat(data, BF_CNOonly_prop_data, 2);
data = padmat(data, BF_Hm4di_prop_data, 2);
[p, tbl, stats] = kruskalwallis(data)
kruskalwallis_Ra = multcompare(stats)

