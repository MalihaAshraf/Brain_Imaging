 %% CC MD 

 close all

files = dir('datasets/short_list_15_Aug_2020/*.xlsx');
% types = {'MD'};

rep = true;             % Generate pdf report true or false
doc = true; 
m_w = true;             % motion weight

for dd = 1:length(files)

close all
data_table = readtable(fullfile(files(dd).folder, files(dd).name));
[~, d_s, ~] =  fileparts(fullfile(files(dd).name));
folder = files(dd).folder;

if m_w
    rep_name = ['figs/', d_s, '_motion_weight_', date]; 
    doc_name = fullfile(folder, [d_s, '_motion_weight_binData.csv']);
else
    rep_name = ['figs/', d_s, '_', date]; 
    doc_name = fullfile(folder, [d_s, '_binData.csv']);
end


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


regions = data_table.Properties.VariableNames(5:end);
subs = unique(data_table.SUBID);
subs_id = data_table.SUBID;
hp_logic = data_table.headphone;

if doc
    bin_data = array2table(zeros(0, length(regions)+1));
    bin_data.Properties.VariableNames = ['bin', regions];
end


mean_per_subj = true;   % Mean of all subjects per region
norm = true;           % Normalize data, true or false
fit_ = true;           % Fit data, true or false
log_ = true;            % No. of days in log, true or false
e_bar = false;           % Error bar, true or false
ind_sub_plot = false;



j = 1;
k = 1;

% for tt = 1:numel(types)
% 	type = types{tt};
    
    if rep
        chap = Chapter([d_s]);
        add(rpt,chap);
        rpt.Document.CurrentPageLayout.PageMargins.Left = '0.5in';
        rpt.Document.CurrentPageLayout.PageMargins.Right = '0.5in';
    end
        
    for rr = 1:numel(regions)
        for hp = 2
            
        if rep
            clf(h)
        else
%             figure('units','normalized','outerposition',[0 0 1 1])
            figure
        end
        
        days_all = [];
		data_all = []; 
        
        for ss = 1:length(subs)
            
            if hp == 1
                days_v = data_table.DAYS(ss == subs_id);
                data = data_table.(regions{rr})(ss == subs_id);
                hp_lab = 'with headphone';
            else 
                days_v = data_table.DAYS(ss == subs_id & ~hp_logic);
                motion = data_table.motion(ss == subs_id & ~hp_logic);
                data = data_table.(regions{rr})(ss == subs_id & ~hp_logic);
            end
            
            days_v(days_v < 0) = 0;
            d_l = 1;    % adjustment factor for log values
            if log_
                days_v = log2(days_v+d_l);
            end
            days_max = log2(200+d_l);
            days_min = log2(-2+d_l);
            
            norm_param = mean(data);
            if norm                
                data = (data)./norm_param; 
            end
            
            if m_w
                w = (1-motion);
                hp_lab = 'without headphone & weighted with motion';
            else
                w = ones(size(data));
                hp_lab = 'without headphone';
            end
            
            if ss == 1
               days_all = days_v;
               data_all = data;
               w_all = w;
            else 
                days_all = [days_all; days_v];
                data_all = [data_all; data];
                w_all = [w_all; w];
                [days_all, ind] = sort(days_all);
                data_all = data_all(ind);
                w_all = w_all(ind);
                clear ind
            end
            
        end
   
        if mean_per_subj
            days_all = days_all(~isnan(data_all));
            data_all = data_all(~isnan(data_all));
            w_all = w_all(~isnan(data_all));
            
            [days_all_u,ia,ic] = unique(ceil(days_all));
            data_all_u = ones(size(days_all_u)).*NaN;
            data_sd_u = ones(size(days_all_u)).*NaN;
            n_u = ones(size(days_all_u)).*NaN;
            
            for uu = 1:length(days_all_u)
                [data_all_u(uu), data_sd_u(uu), ~, ~] = weighted_mean_std(data_all(ic == uu), w_all(ic == uu));
                n_u(uu) = length(find(ic == uu));
            end
            
            if doc
                if rr ==1
                    bin_data(1:length(data_all_u),1) = array2table(days_all_u);
                    bin_data(length(data_all_u)+1:length(data_all_u)*2,1) = array2table(days_all_u);
                    
                    bin_data(1:length(data_all_u),rr+1) = array2table(data_all_u);
                    bin_data(length(data_all_u)+1:length(data_all_u)*2,rr+1) = array2table(data_sd_u);
                else 
                    bin_data(1:length(data_all_u),rr+1) = array2table(data_all_u);
                    bin_data(length(data_all_u)+1:length(data_all_u)*2,rr+1) = array2table(data_sd_u);
                end
            end
            
            if fit_
                ind = 1:6;
				[y_out, mean_y, sign_y, peak_y, min_y] = process_y_bigauss(data_all_u(ind));
                ft = @(a, b, c1, c2, x) bigauss(a, b, c1, c2, x);
                 [f, gof] = fit(days_all_u(ind), y_out, ft,...
                        'StartPoint', [peak_y-min_y, 1.5, 1, 1 ],...
                        'Upper', [Inf, 3, 3, 3],...
                        'Lower', [0, 0, 0.5, 0.5],...
                        'Weights', data_sd_u(ind));   
                    days_f = days_all_u(1):0.1:days_all_u(end);
                    y = feval(f, days_f);
                    y = ((y+min_y).*sign_y)+mean_y;
                    plot(days_f, y', 'Color', 'r', 'LineStyle', '--');
                    hold on
                    errorbar(days_all_u, data_all_u, data_sd_u./2, 'Vertical', '-',...
                        'MarkerSize',5, 'Color', 'k', 'MarkerFaceColor', 'k')
           
            else
%                 scatter(days_all_u, data_all_u, 20, 'k', 'filled')
                errorbar(days_all_u, data_all_u, data_sd_u, 'Vertical', '-',...
                        'MarkerSize',5, 'Color', 'k', 'MarkerFaceColor', 'k')
            end
                        
             if ~ind_sub_plot
    %             xticks = [-5 0 10 30 60 100 200];
                xticks_l = (linspace(days_all_u(1), days_all_u(end), 9));
                xticks = round(2.^(xticks_l)-d_l);
                xlim([days_all_u(1)-0.5, days_all_u(end)+0.5])
                set(gca, 'XTick', xticks_l);
                set(gca', 'XTickLabel', num2str(xticks'));
                xlabel('Scanning days')
                ylabel('Mean diffusivity')
                grid on
             end
        end
        
        suptitle({[strrep(d_s, '_', ' '), ' ', hp_lab], strrep(regions{rr}, '_', ' ')});
       
        
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
    end
    
%     subplot(1, numel(types)+1, numel(types)+1)
%     suptitle(['Region ', num2str(regions(rr))]);
    
% end



% if doc_
%    writetable(T,doc_name) 
%    writetable(T_all, doc2_name)
% end

if rep
    close(rpt);
end

if doc 
    writetable(bin_data, doc_name);
end

end

if rep
close all
h=figure('visible','on');
end
