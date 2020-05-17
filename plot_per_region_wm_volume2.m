%%

types = {'AD', 'FA', 'MD', 'RD'};
n_visits  =  [16 14 14 13 12 13 13 14 13 13 13 11 12 11 10];    % No of visits for each subject
n_sub = length(n_visits);
load days2.mat                       % Day number for each visit for each participant
% days(2, :) = -2;
% regions = [ 4 5 6 23 24 25 26 27 28 29 42 43 44 45];
regions = 1:68;
colors = lines(n_sub);
data_folder = 'datasets/WMV/';

mean_per_subj = false;   % Mean of all subjects per region
norm = true;           % Normalize data, true or false
fit_ = true;           % Fit data, true or false
log_ = true;            % No. of days in log, true or false
e_bar = false;           % Error bar, true or false
ind_sub_plot = false;

labels = readtable('WMV_label_list.xlsx');

rep = true;             % Generate pdf report true or false
rep_name = 'figs/WMV_noMask_may_17th_20';

doc_ = true;
doc_name = 'tables/WMV_noMask_stats.xlsx';
doc2_name = 'tables/WMV_noMask_stats_all.xlsx';

if rep
    import mlreportgen.report.*
    import mlreportgen.dom.*
    rpt = Report(rep_name,'pdf');
    chap = Chapter('Linear fitting');
    add(rpt,chap);
    rpt.Document.CurrentPageLayout.PageMargins.Left = '0.5in';
    rpt.Document.CurrentPageLayout.PageMargins.Right = '0.5in';
end

if rep
    h=figure('visible','off');
%     figure('units','normalized','outerposition',[0.1 0 0.9 1])
	h.PaperUnits = 'inches';
    h.PaperSize = [7.5 4.5];
end

table_vars = [["Regions", "int16"];...
    ["Diff_type", "string"];...
    ["Subject", "int16"];...
    ["n", "int16"];...
    ["RMSE", "double"];...
    ["R_squared", "double"];...
    ["Adj_R_squared", "double"];...
    ["F-statistic", "double"];...
    ["p-value", "double"];...
    ["intercept", "double"];...
    ["slope", "double"]];
table2_vars = [["Regions", "int16"];...
    ["Region_name", "string"];...
    ["Diff_type", "string"];...
    ["T-test", "double"];...
    ["ci_low", "double"];...
    ["ci_high", "double"];...
    ["p-value-test", "double"];...
    ["RMSE", "double"];...
    ["R_squared", "double"];...
    ["Adj_R_squared", "double"];...
    ["F-statistic", "double"];...
    ["p-value-fit", "double"];...
    ["intercept", "double"];...
    ["slope", "double"]];

clear T T_all
T = table('Size',[numel(regions)*numel(types)*n_sub,size(table_vars,1)],... 
	'VariableNames', table_vars(:,1),...
	'VariableTypes', table_vars(:,2));
T_all = table('Size',[numel(regions)*numel(types),size(table2_vars,1)],... 
	'VariableNames', table2_vars(:,1),...
	'VariableTypes', table2_vars(:,2));
j = 1;
k = 1;

for rr = 1:numel(regions)
    r_label =  split(labels.labelName{rr});
    for tt = 1:numel(types)
        clear ds

        type = types{tt};
         
        if rep
            clf(h)
        else
%             figure('units','normalized','outerposition',[0 0 1 1])
            figure
        end
        
		days_all = [];
		data_all = []; 
        
		slopes = [];
        
        for ss = 1:n_sub
            if (~mean_per_subj) && (ind_sub_plot)
                subplot(3, 5, ss);
            end
            
            n_v = n_visits(ss);
            days_v = days(:, ss);
            days_v = days_v(~isnan(days_v));
%            sel = (days_v > 0) & (days_v < 120);
%             days_v = days_v(3:end, :);
%            days_v = days_v(sel, :);
			days_v(days_v < -1) = -1;
            d_l = 2;    % adjustment factor for log values
            if log_
                days_v = log10(days_v+d_l);
            end
            days_max = log10(200+d_l);
            days_min = log10(-2+d_l);
			
            load ([data_folder, 'WMV/', types{tt}, '/S', num2str(ss, '%d'), '_', types{tt}, '.mat']);
            
            if 0%tt == 2
                data = meanFADifValue(rr, 1:length(days_v))';
                if ~sum(data)
                    continue
                end
            else
                data = meanDifValue(rr, 1:length(days_v))';
%                 data_sd = stdDifValue(rr, :);
            end
           
		   norm_param = mean(data);
            d = 0.2;
            if norm                
                data = (data)./norm_param; 
                data_p = data+d*(ss-1);
%                 data_sd = data_sd./norm_param;
            end
            
                if e_bar
                    stdshade(cat(2, data, data_sd),0.2, colors(ss, :), days_v);
    %                 errorbar(days_v, data_p', data_sd, 'Vertical', '.',...
    %                     'MarkerSize',3, 'Color',colors(ss, :), 'MarkerFaceColor',colors(ss, :))
                else
                    scatter(days_v, data', 5, colors(ss, :), 'filled');    
                end
                hold on
                
                if fit_
				
					mdl = fitlm(days_v, data);
                    [p,F] = coefTest(mdl);
                    T(j, :) = {rr, type, ss, mdl.NumObservations,...
                        mdl.RMSE, mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted,...
                        F, p, mdl.Coefficients.Estimate(1), mdl.Coefficients.Estimate(2)};
                    slopes = [slopes;  mdl.Coefficients.Estimate(2)];

                    plot(days_v, mdl.Fitted, 'Color', colors(ss, :), 'LineStyle', '--');
                else
				   hh(ss) = plot(days_v, data', 'Color', colors(ss, :));
                end
                hold on
				j = j+1;

    %             label{ss} = ['Subject ' num2str(ss)];
    %            title(['S' num2str(ss)])

            
%                 hh(ss) = scatter(days_v, data', 5, colors(ss, :), 'filled');  
                hold on
                if ss == 1
                   days_all = days_v;
                   data_all = data;
                else 
                    days_all = [days_all; days_v];
                    data_all = [data_all; data];
                    [days_all, ind] = sort(days_all);
                    data_all = data_all(ind);
                    clear ind
                end
              
        end
        
%        if mean_per_subj
            [days_all_u,ia,ic] = unique(days_all);
            data_all_u = ones(size(days_all_u)).*NaN;
            data_sd_u = ones(size(days_all_u)).*NaN;
            n_u = ones(size(days_all_u)).*NaN;
            
            for uu = 1:length(days_all_u)
                data_all_u(uu) = median(data_all(ic == uu));
                data_sd_u(uu) = iqr(data_all(ic == uu));
                n_u(uu) = length(find(ic == uu));
            end
            
            if fit_
				mdl = fitlm(days_all_u, data_all_u);
                [p,F] = coefTest(mdl);
                [h1,p2,ci,stats] = ttest(slopes);
                T_all(k, :) = {rr, r_label{2}, type,  h1, ci(1), ci(2), p2,...
                    mdl.RMSE, mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted,...
                    F, p, mdl.Coefficients.Estimate(1), mdl.Coefficients.Estimate(2)};
                k = k+1;

               plot(days_all_u, mdl.Fitted,... 
                       'LineWidth', 2, 'LineStyle', '-', 'Color', [0 0 0]);
					   
                %text(days_all_u(end-15), max(data_all)*1.2,... 
                 %      ['rmse: ', num2str(gof.rmse)]);

             end
                        
             
%             xticks = [-5 0 10 30 60 100 200];
            xticks_l = (linspace(days_all_u(1), days_all_u(end), 10));
            xticks = round(10.^(xticks_l)-d_l);
            set(gca, 'XTick', xticks_l);
            set(gca', 'XTickLabel', num2str(xticks'));
            xlabel('Scanning days')
            ylabel('Mean diffusivity')
            grid on
%        end
        
        suptitle([type, ' ', r_label{2}]);
        
%         ylim([0.8 d*20])

%          if tt == 4
%             legend(hh, label, 'location', 'eastoutside')
%             legend boxoff
%         end
    
        if rep
%            fig = Figure(gcf);
%            fig.Width = '8in';
%            fig.Height = '5in';
%            fig.Scaling = 'custom';
            add(rpt, Figure(h));
        end
    
    end
    
%     subplot(1, numel(types)+1, numel(types)+1)
%     suptitle(['Region ', num2str(regions(rr))]);
    
end

if rep
    close(rpt);
end

if doc_
   writetable(T,doc_name) 
   writetable(T_all, doc2_name)
end

close all
h=figure('visible','on');