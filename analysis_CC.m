%% Whole brain

types = {'FA', 'MD', 'RD', 'L1', 'L2', 'L3'};
n_visits  =  [16 14 14 13 12 13 13 14 13 13 13 11 12 11 10];
n_sub = length(n_visits);
load days2.mat
load b_data.mat
% days(2, :) = -2;
% regions = [ 4 5 6 23 24 25 26 27 28 29 42 43 44 45];
regions = 1:5;
colors = lines(n_sub);
rep = false;
log_ = true;

%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw DATA

norm = true;
fit_ = false;
path = 'CC_DTIdata\CC_DTIdata\';

if rep
    import mlreportgen.report.*
    import mlreportgen.dom.*
    rpt = Report('CC_w_behavioral','pdf');
    chap = Chapter('Raw data');
    add(rpt,chap);
    rpt.Document.CurrentPageLayout.PageMargins.Left = '0.5in';
    rpt.Document.CurrentPageLayout.PageMargins.Right = '0.5in';
end

% figure('units','normalized','outerposition',[0 0 1 1]);

for tt = 1:numel(types)
    figure
    
    type = types{tt};
    for rr = 1:numel(regions)

        i = 1;
        subplot(2, 3, rr)
        axis square
        
        for ss = 1:n_sub
            
            fileID = [path, 'CC_' type, 'files\S', num2str(ss, '%0.2d'), type, '_CC.mat'];
            load (fileID);
            
               n_v = n_visits(ss);
               eval(['data = meanDifValue(rr, :);']);
               norm_param = mean(data);
               
               if norm
                  data = (data)./norm_param; 
               end
               
               days_v = days(:, ss);
               days_v = days_v(~isnan(days_v));
%                days_v = days_v(2:end, :);
               if log_
                   days_v = log10(days_v+50);
               end
               
               if length(days_v) ~= length(data)
                    break
                end
               scatter(days_v, data', 5, colors(ss, :), 'filled')
               hold on
               
               if fit_
                   ft = fittype('(a-c)/(1+exp(-b*(x-d)))+c');
%                    f = fit((days_v), data, ft, 'Weights', data_sd,...
%                        'Upper', [Inf, 0, Inf Inf]);
                   f = fit((days_v), data, ft, 'Weights', data_sd);
                   y = feval(f, days_v(1):1:days_v(end));
                   hh(rr) = plot(days_v(1):1:days_v(end), y', 'Color', colors(ss, :));
               else
                   hh(rr) = plot(days_v, data', 'Color', colors(ss, :));
               end
               
               label{rr} = ['Region ', num2str(regions(rr))];
               hold on
               i = i+ n_v;
        end
        xlim([1.5 2.5])
        title(['Region ' num2str(regions(rr))])

        %         legend(hh, label);
%         legend boxoff
        
    end 
    
    if norm
        if fit_
            suptitle(['Fitted normalized ' type])
        else
            suptitle(['Normalized ' type])
        end  
    else
        suptitle(type)
    end
    
    if rep
        fig = Figure(gcf);
        fig.Width = '8in';
        fig.Height = '5in';
        fig.Scaling = 'custom';
        add(rpt, fig);
%         delete(gcf);
    end
end
