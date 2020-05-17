%%

types = {'AD', 'FA', 'MD', 'RD'};
n_visits  =  [16 14 14 13 12 13 13 14 13 13 13 11 12 11 10];    % No of visits for each subject
n_sub = length(n_visits);
load days.mat                       % Day number for each visit for each participant
% days(2, :) = -2;
% regions = [ 4 5 6 23 24 25 26 27 28 29 42 43 44 45];
regions = 1:48;
colors = lines(n_sub);
data_folder = 'datasets/DTI/';

mean_all_subj = false;   % Mean of all subjects per region
norm = true;           % Normalize data, true or false
fit_ = true;           % Fit data, true or false
log_ = true;            % No. of days in log, true or false
e_bar = false;           % Error bar, true or false

rep = true;             % Generate pdf report true or false
rep_name = 'figs/DTI_analysis_may_17th_20';

doc_ = true;
doc_name = 'tables/DTI_stats.xlsx';
doc2_name = 'tables/DTI_stats_all.xlsx';

if rep
    import mlreportgen.report.*
    import mlreportgen.dom.*
    rpt = Report(rep_name,'pdf');
    chap = Chapter('Linear fitting per subject');
    add(rpt,chap);
    rpt.Document.CurrentPageLayout.PageMargins.Left = '0.5in';
    rpt.Document.CurrentPageLayout.PageMargins.Right = '0.5in';
end

if rep
   h=figure('visible','off');
%     set(h,'visible','off');
%     figure_handle = h;
%     figure_handle.WindowState = 'minimized';
    h.PaperUnits = 'inches';
    h.PaperSize = [7.5 4.5];
%     h.Width = '8in';
%     h.Height = '5in';
%     h.Scaling = 'custom';

%     figure('units','normalized','outerposition',[0.1 0 0.9 1])
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
T = table('Size',[numel(regions)*numel(types)*n_sub,size(table_vars,1)],... 
	'VariableNames', table_vars(:,1),...
	'VariableTypes', table_vars(:,2));
T_all = table('Size',[numel(regions)*numel(types),size(table2_vars,1)],... 
	'VariableNames', table2_vars(:,1),...
	'VariableTypes', table2_vars(:,2));
j = 1;
k = 1;

for rr = 1:numel(regions)

    for tt = 1:numel(types)
        clear ds

        type = types{tt};
         
         i = 1;
        if rep
%             set(0, 'CurrentFigure', figure_handle);
            clf(h)
        else
%             figure('units','normalized','outerposition',[0 0 1 1])
            figure
        end
        
       days_all = [];
	   data_all = []; 
        
	   slopes = [];
        for ss = 1:n_sub
            if (~mean_all_subj) && (ind_sub_plot)
                subplot(3, 5, ss);
            end
            
            n_v = n_visits(ss);
            days_v = days(:, ss);
            days_v = days_v(~isnan(days_v));
            days_v(days_v < -1) = -1;
%             days_v = days_v(2:end, :);
            d_l = 2;    % adjustment factor for log values
            if log_
                days_v = log10(days_v+d_l);
            end
            days_max = log10(200+d_l);
            days_min = log10(-2+d_l);
    
            fileID = fopen([data_folder, type '/' type '/ske', num2str(regions(rr)) '_' type '_values.txt'], 'r');
            formatSpec = '%f %f %f %f';
            sz = [4 Inf];
            ds = fscanf(fileID, formatSpec, sz);
            ds = ds';
            fclose(fileID);
            
            data = ds(i:i+n_v-1, 1);
            data_sd = ds(i:i+n_v-1, 2);
            norm_param = mean(data);
            
            d = 0.2;
            if norm
                data = (data)./norm_param; 
                data_p = data+d*(ss-1);
                data_sd = data_sd./norm_param;
            end
            
                if e_bar
                    stdshade(cat(2, data, data_sd),0.2, colors(ss, :), days_v);
                else
                    scatter(days_v, data', 5, colors(ss, :), 'filled');    
                end
                hold on
                
                if fit_

                    mdl = fitlm(days_v, data, 'Weights', data_sd);
                    [p,F] = coefTest(mdl);
                    T(j, :) = {rr, type, ss, mdl.NumObservations,...
                        mdl.RMSE, mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted,...
                        F, p, mdl.Coefficients.Estimate(1), mdl.Coefficients.Estimate(2)};
                    slopes = [slopes;  mdl.Coefficients.Estimate(2)];

                       plot(days_v, mdl.Fitted, 'Color', colors(ss, :), 'LineStyle', '--');
                else
                       hh(ss) = plot(days_v, data', 'Color', colors(ss, :));
                end

                i = i+ n_v;
                j = j+1;
                
                if ind_sub_plot
                    title(['S' num2str(ss)])
                    xticks_l = (linspace(days_v(1), days_v(end), 4));
                    xticks = round(10.^(xticks_l)-d_l);
                    set(gca, 'XTick', xticks_l);
                    set(gca', 'XTickLabel', num2str(xticks'));
                    grid on
                end
                
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
        
%         if mean_all_subj
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
                mdl = fitlm(days_all_u, data_all_u, 'Weights', data_sd_u);
                [p,F] = coefTest(mdl);
                [h1,p2,ci,stats] = ttest(slopes);
                T_all(k, :) = {rr, type,  h1, ci(1), ci(2), p2,...
                    mdl.RMSE, mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted,...
                    F, p, mdl.Coefficients.Estimate(1), mdl.Coefficients.Estimate(2)};
                k = k+1;

               plot(days_all_u, mdl.Fitted,... 
                       'LineWidth', 2, 'LineStyle', '-', 'Color', [0 0 0]);

             end
            
            xticks_l = (linspace(days_all_u(1), days_all_u(end), 10));
            xticks = round(10.^(xticks_l)-d_l);
            set(gca, 'XTick', xticks_l);
            set(gca', 'XTickLabel', num2str(xticks'));
            xlabel('Scanning days')
            ylabel('Mean diffusivity')
            grid on
        
        suptitle([type, ' Region ', num2str(regions(rr))]);
    
        if rep
%             fig = Figure(h);
%             fig.Width = '8in';
%             fig.Height = '5in';
%             fig.Scaling = 'custom';
            add(rpt, Figure(h));
        end
    
    end
        
end

if rep
    close(rpt);
end

if doc_
   writetable(T,doc_name) 
   writetable(T_all, doc2_name)
end