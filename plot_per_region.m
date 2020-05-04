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
ind_sub_plot = false;   % subplot for individual subjects
norm = true;           % Normalize data, true or false
fit_ = true;           % Fit data, true or false
log_ = true;            % No. of days in log, true or false
e_bar = false;           % Error bar, true or false

rep = true;             % Generate pdf report true or false
rep_name = 'figs/DTI_analysis_may_4th_20';

doc_ = true;
doc_name = 'tables/DTI_stats.xlsx';

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

T = table('Size',[numel(regions)*numel(types)*n_sub,size(table_vars,1)],... 
	'VariableNames', table_vars(:,1),...
	'VariableTypes', table_vars(:,2));
j = 1;

for rr = 1:numel(regions)

    for tt = 1:numel(types)
        clear ds
%         subplot(1, n_sub, ss+(tt-1)*n_sub)
%         subplot(1, 4, tt)

        type = types{tt};
         
         i = 1;
        if rep
%             set(0, 'CurrentFigure', figure_handle);
            clf(h)
        else
%             figure('units','normalized','outerposition',[0 0 1 1])
            figure
        end
        
        if mean_all_subj
           days_all = [];
           data_all = []; 
        end
        
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
                data = (data - norm_param)./norm_param; 
                data_p = data+d*(ss-1);
                data_sd = data_sd./norm_param;
            end
            
            if ~mean_all_subj
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
                    
    %                    f = fit(days_v, data, 'gauss2', 'Weights', data_sd);
%                        y = feval(f, days_v(1):0.001:days_v(end));
%                        y_p = y+d*(ss-1);
                       plot(days_v, mdl.Fitted, 'Color', colors(ss, :), 'LineStyle', '--');
                else
                       hh(ss) = plot(days_v, data', 'Color', colors(ss, :));
                end
%                 hold on

    %             label{ss} = ['Subject ' num2str(ss)];
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
                
               
            else
                hh(ss) = scatter(days_v, data', 5, colors(ss, :), 'filled');  
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
        end
        
        if mean_all_subj
            [days_all_u,ia,ic] = unique(days_all);
            data_all_u = ones(size(days_all_u)).*NaN;
            data_sd_u = ones(size(days_all_u)).*NaN;
            n_u = ones(size(days_all_u)).*NaN;
            
            for uu = 1:length(days_all_u)
                data_all_u(uu) = median(data_all(ic == uu));
                data_sd_u(uu) = iqr(data_all(ic == uu));
                n_u(uu) = length(find(ic == uu));
            end
            
            if 1 % choose between error region or error bar
                stdshade(cat(2, data_all_u, data_sd_u),0.2, [0.5, 0.5, 0.5], days_all_u);
%                 hold on
%                 scatter(days_all, data_all, 3, 'filled');
            else
                errorbar(days_all_u, data_all_u, data_sd_u);
            end
            
            if fit_
                ft = '(p*(a1*exp(-((x-b1)/c1)^2)+d1)) + ((1-p)*(a2*exp(-((x-b2)/c2)^2)+d2))';
                [~, ind] = max(abs(data_all_u - mean(data_all_u)));
                if data_all_u(ind) < 1
                    st_pts = [-0.5, -0.5, 1.1, 2.1, 0.1, 0.5, 0.5, 1, 0.5];
                else
                    st_pts = [0.5, 0.5, 1.1, 2.1, 0.1, 0.5, 0.5, 1, 0.5];
                end

                f = fit(days_all_u, data_all_u, ft,...
                    'StartPoint', st_pts);

%                    f = fit(days_v, data, 'gauss2', 'Weights', data_sd);
                   y = feval(f, days_all_u(1):0.001:days_all_u(end));
                   y_p = y+d*(ss-1);
                   hh(ss) = plot(days_all_u(1):0.001:days_all_u(end), y',... 
                       'LineStyle', '--', 'Color', [0.5, 0.5, 0.5]);

             end
            
%             xticks = [-5 0 10 30 60 100 200];
            xticks_l = (linspace(days_all_u(1), days_all_u(end), 10));
            xticks = round(10.^(xticks_l)-10);
            set(gca, 'XTick', xticks_l);
            set(gca', 'XTickLabel', num2str(xticks'));
            xlabel('Scanning days')
            ylabel('Median z-score normalized data (error region: iqr)')
            grid on
        else
             if ~ind_sub_plot
                    xlim([days_min days_max]);
                    xticks_l = (linspace(0, days_max, 10));
                    xticks = round(10.^(xticks_l)-d_l);
                    set(gca, 'XTick', xticks_l);
                    set(gca', 'XTickLabel', num2str(xticks'));
                    xlabel('Scanning days')
                    ylabel('Mean diffusivity')
                    grid on
             end
        end
        
        suptitle([type, ' Region ', num2str(regions(rr))]);
    
        if rep
%             fig = Figure(h);
%             fig.Width = '8in';
%             fig.Height = '5in';
%             fig.Scaling = 'custom';
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
end