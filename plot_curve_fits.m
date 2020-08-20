 %% ASEG 

 close all

d_s = 'ASEG';

types = {'FA', 'MD', 'RD', 'AD'};
% types = {'FA'};
n_visits  =  [16 14 14 13 12 13 13 14 13 13 13 11 12 11 10];    % No of visits for each subject
n_sub = length(n_visits);
load days2.mat                       % Day number for each visit for each participant
% days(2, :) = -2;
% regions = [ 4 5 6 23 24 25 26 27 28 29 42 43 44 45];
% regions = [8, 13, 14, 27, 28, 29, 30, 31];
colors = lines(n_sub);

switch d_s
    case 'ASEG'
        data_folder = 'datasets/ASEG/';
        labels = readtable('ASEG_label_list.xlsx');
        regions = [37, 38, 39, 40, 41];

    case 'CC'
        data_folder = 'datasets_old/CC_DTIdata/';
        regions = 1:5;
end

mask_label = {' (No mask)', ' (FA mask)'};

mean_per_subj = true;   % Mean of all subjects per region
norm = true;           % Normalize data, true or false
fit_ = true;           % Fit data, true or false
log_ = true;            % No. of days in log, true or false
e_bar = false;           % Error bar, true or false
ind_sub_plot = false;

rep = true;             % Generate pdf report true or false
rep_name = ['figs/', d_s, '_fits_', date]; 


% doc_ = true;
% doc_name = 'tables/ASEG_noMask_stats.xlsx';
% doc2_name = 'tables/ASEG_noMask_stats_all.xlsx';

if rep
    import mlreportgen.report.*
    import mlreportgen.dom.*
    rpt = Report(rep_name,'pdf');
    
end

if rep
    h=figure('visible','off');
    h.PaperUnits = 'inches';
    h.PaperSize = [7.5 4.5];
%     figure('units','normalized','outerposition',[0.1 0 0.9 1])
end


j = 1;
k = 1;

for ff = 1:2
for tt = 1:numel(types)
	type = types{tt};
    
    if rep
        chap = Chapter([types{tt}, mask_label{ff}]);
        add(rpt,chap);
        rpt.Document.CurrentPageLayout.PageMargins.Left = '0.5in';
        rpt.Document.CurrentPageLayout.PageMargins.Right = '0.5in';
    end
        
    for rr = 1:numel(regions)
        clear ds
        
        if strcmp(d_s, 'ASEG')
            r_label =  labels.labelName{regions(rr)};
        end
        
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
            
            days_v = days(:, ss);
            days_v = days_v(~isnan(days_v));
%            sel = (days_v >= 0);% & (days_v < 120);
%             days_v = days_v(3:end, :);
            days_v(days_v < 0) = 0;
            d_l = 1;    % adjustment factor for log values
            if log_
                days_v = log2(days_v+d_l);
            end
            days_max = log2(200+d_l);
            days_min = log2(-2+d_l);
			
            switch d_s
                case 'ASEG'            
                    load ([data_folder, 'ASEG/', types{tt}, '/S', num2str(ss, '%d'), '_', types{tt}, '.mat']);
                case 'CC'
                    load ([data_folder, 'CC_', types{tt}, 'files/S', num2str(ss, '%.2d'), types{tt}, '_CC.mat']);
            end
            if ff == 2%tt == 2
                data = meanFADifValue(regions(rr), 1:length(days_v))';
                if ~sum(data)
                    continue
                end
                if strcmp(d_s, 'CC')
                    data_sd = [];
                else
                    data_sd = stdFADifValue(regions(rr), 1:length(days_v));
                end 
            else
                data = meanDifValue(regions(rr), 1:length(days_v))';
                if ~sum(data)
                    continue
                end
                if strcmp(d_s, 'CC')
                    data_sd = [];
                else
                    data_sd = stdDifValue(regions(rr), 1:length(days_v));
                end 
            end
           
            if 0
               days_v = days_v (1:10);
               data = data(1:10);
               data_sd = data_sd(1:10);
            end
            
			norm_param = mean(data);
            d = 0;
            if norm                
                data = (data)./norm_param; 
                data_p = data+d*(ss-1);
                data_sd = data_sd./norm_param;
            end
            
            
            [y_out, mean_y, sign_y, peak_y, min_y] = process_y_bigauss(data);
            if ~mean_per_subj
                if e_bar
%                     stdshade(cat(2, data, data_sd),0.2, colors(ss, :), days_v);
                    errorbar(days_v, data_p', data_sd, 'Vertical', '.',...
                        'MarkerSize',3, 'Color',colors(ss, :), 'MarkerFaceColor',colors(ss, :))
                else
                    
                    scatter(days_v, data, 20, colors(ss, :), 'filled');    
                end
                hold on
                
                if fit_
% 					mdl = fitlm(days_v, data);
%                     [p,F] = coefTest(mdl);
%                     T(j, :) = {rr, type, ss, mdl.NumObservations,...
%                         mdl.RMSE, mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted,...
%                         F, p, mdl.Coefficients.Estimate(1), mdl.Coefficients.Estimate(2)};
%                     slopes = [slopes;  mdl.Coefficients.Estimate(2)];
                    
%                     [pol,S] = polyfit(days_v, data, 3);
%                     plot(days_v, polyval(pol, days_v), 'Color', colors(ss, :), 'LineStyle', '--');
                    
                    [y_out, mean_y, sign_y, peak_y, min_y] = process_y_bigauss(data);
                    ft = @(a, b, c1, c2, x) bigauss(a, b, c1, c2, x);
                    
%                     f = fit(days_v, y_out, ft);
                    
                    [f, gof] = fit(days_v, y_out, ft,...
                        'StartPoint', [peak_y-min_y, 1.2, 1, 1 ],...
                        'Upper', [Inf, 2, 2, 2],...
                        'Lower', [0, 0, 0.5, 0.5],...
                        'Weights', data_sd);
                    
                    days_f = days_v;
%                     days_f = days_v(1):0.001:days_v(end);

                    y = feval(f, days_f);
                    y = ((y+min_y).*sign_y)+mean_y;
                    
%                     if sign_y < 0
%                        2 
%                     end
                    hh(ss) = plot(days_f, y', 'Color', colors(ss, :));
                    
                else
                       hh(ss) = plot(days_v, y_out', 'Color', colors(ss, :));
                end
                hold on
				j = j+1;

    %             label{ss} = ['Subject ' num2str(ss)];
    %            title(['S' num2str(ss)])

            else
            %    hh(ss) = scatter(days_v, data', 5, colors(ss, :), 'filled');  
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
            
            if ind_sub_plot
%             xticks = [-5 0 10 30 60 100 200];
                xticks_l = (linspace(days_v(1), days_v(end), 5));
                xticks = round(2.^(xticks_l)-d_l);
                set(gca, 'XTick', xticks_l);
                set(gca', 'XTickLabel', num2str(xticks'));
                if ss == 13
                    xlabel('Scanning days')
                end
                if ss == 6
                    ylabel('Mean diffusivity')
                end
                title(['rmse: ', num2str(gof.rmse)])
                grid on
             end
        end
        
        if mean_per_subj
            [days_all_u,ia,ic] = unique(ceil(days_all));
            data_all_u = ones(size(days_all_u)).*NaN;
            data_sd_u = ones(size(days_all_u)).*NaN;
            n_u = ones(size(days_all_u)).*NaN;
            
            for uu = 1:length(days_all_u)
                data_all_u(uu) = mean(data_all(ic == uu));
                data_sd_u(uu) = std(data_all(ic == uu))./sqrt(length(data_all(ic == uu)));
                n_u(uu) = length(find(ic == uu));
            end
            
            if fit_
                if 0
				mdl = fitlm(days_all_u, data_all_u);
                [p,F] = coefTest(mdl);
                [h1,p2,ci,stats] = ttest(slopes);
                T_all(k, :) = {rr, r_label, type,  h1, ci(1), ci(2), p2,...
                    mdl.RMSE, mdl.Rsquared.Ordinary, mdl.Rsquared.Adjusted,...
                    F, p, mdl.Coefficients.Estimate(1), mdl.Coefficients.Estimate(2)};
                k = k+1;

               plot(days_all_u, mdl.Fitted,... 
                       'LineWidth', 2, 'LineStyle', '-', 'Color', [0 0 0]);
					   
                %text(days_all_u(end-15), max(data_all)*1.2,... 
                 %      ['rmse: ', num2str(gof.rmse)]);
                else
                    ind = 1:6;
                    [y_out, mean_y, sign_y, peak_y, min_y] = process_y_bigauss(data_all_u(ind));
                    ft = @(a, b, c1, c2, x) bigauss(a, b, c1, c2, x);
                     [f, gof] = fit(days_all_u(ind), y_out, ft,...
                            'StartPoint', [peak_y-min_y, 1.5, 1, 1 ],...
                            'Upper', [Inf, 4, 5, 5],...
                            'Lower', [0, 0, 0.5, 0.5],...
                            'Weights', data_sd_u(ind));   
                        days_f = days_all_u(1):0.1:days_all_u(end);
                        y = feval(f, days_f);
                        y = ((y+min_y).*sign_y)+mean_y;
                        plot(days_f, y', 'Color', 'r', 'LineStyle', '--');
                        hold on
                        errorbar(days_all_u, data_all_u, data_sd_u, 'Vertical', '-',...
                            'MarkerSize',5, 'Color', 'k', 'MarkerFaceColor', 'k')

                end
            else
%                 scatter(days_all_u, data_all_u, 20, 'k', 'filled')
                errorbar(days_all_u, data_all_u, data_sd_u, 'Vertical', '-',...
                        'MarkerSize',5, 'Color',colors(ss, :), 'MarkerFaceColor',colors(ss, :))
            end
                        
             if ~ind_sub_plot
    %             xticks = [-5 0 10 30 60 100 200];
                xticks_l = (linspace(days_all_u(1), days_all_u(end), 9));
                xticks = round(2.^(xticks_l)-d_l);
                set(gca, 'XTick', xticks_l);
                set(gca', 'XTickLabel', num2str(xticks'));
                xlabel('Scanning days')
                ylabel('Mean diffusivity')
                grid on
             end
        end
        
        if strcmp(d_s, 'CC')
            suptitle([type, ' Region CC ', num2str(regions(rr))]);
        else
            suptitle([type, ' Region ', num2str(regions(rr)), ' ', labels.labelName{regions(rr)}]);
        end
        
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
end

if rep
    close(rpt);
end

% if doc_
%    writetable(T,doc_name) 
%    writetable(T_all, doc2_name)
% end

if rep
close all
h=figure('visible','on');
end
