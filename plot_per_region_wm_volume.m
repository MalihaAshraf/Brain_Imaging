%%

types = {'AD', 'FA', 'MD', 'RD'};
n_visits  =  [16 14 14 13 12 13 13 14 13 13 13 11 12 11 10];    % No of visits for each subject
n_sub = length(n_visits);
load days2.mat                       % Day number for each visit for each participant
% days(2, :) = -2;
% regions = [ 4 5 6 23 24 25 26 27 28 29 42 43 44 45];
regions = 1:191;
colors = lines(n_sub);

mean_per_subj = true;   % Mean of all subjects per region
norm = true;           % Normalize data, true or false
fit_ = false;           % Fit data, true or false
log_ = true;            % No. of days in log, true or false
e_bar = true;           % Error bar, true or false

rep = true;             % Generate pdf report true or false
rep_name = 'Per_subject_znorm_median_wm_volume';

if rep
    import mlreportgen.report.*
    import mlreportgen.dom.*
    rpt = Report(rep_name,'pdf');
    chap = Chapter('Median z-score normalized data');
    add(rpt,chap);
    rpt.Document.CurrentPageLayout.PageMargins.Left = '0.5in';
    rpt.Document.CurrentPageLayout.PageMargins.Right = '0.5in';
end

if rep
    figure
%     figure('units','normalized','outerposition',[0.1 0 0.9 1])
end

for rr = 1:numel(regions)

    for tt = 1:numel(types)
        clear ds
%         subplot(1, n_sub, ss+(tt-1)*n_sub)
%         subplot(1, 4, tt)

        type = types{tt};
         
         i = 1;
        if rep
            clf
        else
%             figure('units','normalized','outerposition',[0 0 1 1])
            clf
        end
        
        if mean_per_subj
           days_all = [];
           data_all = []; 
        end
        
        for ss = 1:n_sub
            if ~mean_per_subj
                subplot(3, 5, ss);
            end
            
            n_v = n_visits(ss);
            days_v = days(:, ss);
            days_v = days_v(~isnan(days_v));
            days_v = days_v(3:end, :);
            if log_
                days_v = log10(days_v+10);
            end
    
%             fileID = fopen([type '/' type '/ske', num2str(regions(rr)) '_' type '_values.txt'], 'r');
%             formatSpec = '%f %f %f %f';
%             sz = [4 Inf];
%             ds = fscanf(fileID, formatSpec, sz);
%             ds = ds';
%             fclose(fileID);
%             
%             data = ds(i+1:i+n_v-1, 1);
%             data_sd = ds(i+1:i+n_v-1, 2);
%             norm_param = mean(data);
%             
            load (['freesurferDATA/', types{tt}, '/S', num2str(ss, '%.2d'), types{tt}, '_dwiECI.mat']);
            
            if tt == 2
                data = meanFADifValue(rr, 3:end)';
            else
                data = meanDifValue(rr, 3:end)';
%                 data_sd = stdDifValue(rr, :);
            end
           
            d = 0.2;
            if norm                
                data = zscore(data);
                data_p = data+d*(ss-1);
            end
            
            if ~mean_per_subj
                if e_bar
                    stdshade(cat(2, data, data_sd),0.2, colors(ss, :), days_v);
    %                 errorbar(days_v, data_p', data_sd, 'Vertical', '.',...
    %                     'MarkerSize',3, 'Color',colors(ss, :), 'MarkerFaceColor',colors(ss, :))
                else
                    scatter(days_v, data_p', 5, colors(ss, :), 'filled');    
                end
                hold on
                
                if fit_
    %                 ft = '(p*a1*exp(-((x-b1)/c1)^2)+d1) + ((1-p)*a2*exp(-((x-b2)/c2)^2)+d2)';
                    ft = '(p*(a1*exp(-((x-b1)/c1)^2)+d1)) + ((1-p)*(a2*exp(-((x-b2)/c2)^2)+d2))';
                    [~, ind] = max(abs(data - mean(data)));
                    if data(ind) < 1
                        st_pts = [-0.5, -0.5, 1.1, 2.1, 0.1, 0.5, 0.5, 1, 0.5];
                    else
                        st_pts = [0.5, 0.5, 1.1, 2.1, 0.1, 0.5, 0.5, 1, 0.5];
                    end
    %                 f = fit(days_v, data, ft);

                    f = fit(days_v, data, ft,...
                        'StartPoint', st_pts);

    %                    f = fit(days_v, data, 'gauss2', 'Weights', data_sd);
                       y = feval(f, days_v(1):0.001:days_v(end));
                       y_p = y+d*(ss-1);
                       hh(ss) = plot(days_v(1):0.001:days_v(end), y', 'Color', colors(ss, :));
                else
                       hh(ss) = plot(days_v, data', 'Color', colors(ss, :));
                end
                hold on

    %             label{ss} = ['Subject ' num2str(ss)];
                i = i+ n_v;
                title(['S' num2str(ss)])

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
        
        if mean_per_subj
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
            
%             xticks = [-5 0 10 30 60 100 200];
            xticks_l = (linspace(days_all_u(1), days_all_u(end), 10));
            xticks = round(10.^(xticks_l)-10);
            set(gca, 'XTick', xticks_l);
            set(gca', 'XTickLabel', num2str(xticks'));
            xlabel('Scanning days')
            ylabel('Median z-score normalized data (error region: iqr)')
            grid on
        end
        
        suptitle([type, ' Region ', num2str(regions(rr))]);
         
%         ylim([0.8 d*20])

%          if tt == 4
%             legend(hh, label, 'location', 'eastoutside')
%             legend boxoff
%         end
    
        if rep
            fig = Figure(gcf);
            fig.Width = '8in';
            fig.Height = '5in';
            fig.Scaling = 'custom';
            add(rpt, fig);
        end
    
    end
    
%     subplot(1, numel(types)+1, numel(types)+1)
%     suptitle(['Region ', num2str(regions(rr))]);
    
end

if rep
    close(rpt);
end