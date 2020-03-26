
%%

types = {'AD', 'FA', 'MD', 'RD'};
n_visits  =  [16 14 14 13 12 13 13 14 13 13 13 11 12 11 10];
n_sub = length(n_visits);
load days.mat
days(2, :) = -2;
regions = [ 4 5 6 23 24 25 26 27 28 29 42 43 44 45];
colors = lines(numel(regions));
norm = true;
fit_ = true;

for tt = 1:numel(types)
    figure('units','normalized','outerposition',[0 0 1 1])
    type = types{tt};
     i = 1;
     
    for ss = 1:n_sub
        clear ds
%         subplot(1, n_sub, ss+(tt-1)*n_sub)
        subplot(1, n_sub, ss)
        n_v = n_visits(ss);
        
        days_v = days(:, ss);
        days_v = days_v(~isnan(days_v));
        days_v = days_v(2:end, :);
        
        for rr = 1:numel(regions)
            fileID = fopen([type '/' type '/ske', num2str(regions(rr)) '_' type '_values.txt'], 'r');
            formatSpec = '%f %f %f %f';
            sz = [4 Inf];
            ds = fscanf(fileID, formatSpec, sz);
            ds = ds';
            fclose(fileID);
            
            data = ds(i+1:i+n_v-1, 1);
            data_sd = ds(i+1:i+n_v-1, 2);
            norm_param = data(2);
            
            if norm
                data = (data-norm_param)./norm_param; 
                data = data+0.1*(rr-1);
            end
            
            scatter(days_v, data', 5, colors(rr, :), 'filled')
            hold on
            
            if fit_
                    ft = fittype('a/(1+exp(-b*x))');
                   f = fit((days_v), data, ft, 'Weights', data_sd,...
                    'Lower', [-Inf -Inf]);
                   y = feval(f, days_v(1):1:days_v(end));
                   hh(rr) = plot(days_v(1):1:days_v(end), y', 'Color', colors(rr, :));
            else
                   plot(days_v, data', 'Color', colors(rr, :))
            end
            hold on
            
            if tt == 1
               title(['Sub ' num2str(ss)]) 
            end
            if ss == 1
               ylabel(type) 
            end
            label{rr} = ['Region ', num2str(regions(rr))];
        end
        i = i+ n_v;
        
        ylim([-0.2 1.4])

    end
    legend(hh, label)
    legend boxoff
    
end

%%

types = {'AD', 'FA', 'MD', 'RD'};
n_visits  =  [16 14 14 13 12 13 13 14 13 13 13 11 12 11 10];
n_sub = length(n_visits);
load days.mat
% days(2, :) = -2;
regions = [ 4 5 6 23 24 25 26 27 28 29 42 43 44 45];
colors = lines(numel(regions));
norm = true;
fit_ = true;
log_ = true;

rep = false;

if rep
    import mlreportgen.report.*
    import mlreportgen.dom.*
    rpt = Report('Per_participant','pdf');
    chap = Chapter('Fitted data');
    add(rpt,chap);
    rpt.Document.CurrentPageLayout.PageMargins.Left = '0.5in';
    rpt.Document.CurrentPageLayout.PageMargins.Right = '0.5in';
end

i = 1;

for ss = 1:n_sub
%     clf
    figure('units','normalized','outerposition',[0 0 1 1])

    n_v = n_visits(ss);
    days_v = days(:, ss);
    days_v = days_v(~isnan(days_v));
    days_v = days_v(2:end, :);
    if log_
        days_v = log10(days_v+10);
    end
    
    for tt = 1:numel(types)
        clear ds
%         subplot(1, n_sub, ss+(tt-1)*n_sub)
        subplot(1, numel(types), tt)
         type = types{tt};
         
        for rr = 1:numel(regions)
            fileID = fopen([type '/' type '/ske', num2str(regions(rr)) '_' type '_values.txt'], 'r');
            formatSpec = '%f %f %f %f';
            sz = [4 Inf];
            ds = fscanf(fileID, formatSpec, sz);
            ds = ds';
            fclose(fileID);
            
            data = ds(i+1:i+n_v-1, 1);
            data_sd = ds(i+1:i+n_v-1, 2);
            norm_param = mean(data);
            
            if norm
                data = (data)./norm_param; 
                data_p = data+0.1*(rr-1);
            end
            
            scatter(days_v, data_p', 5, colors(rr, :), 'filled')
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
                   y_p = y+0.1*(rr-1);
                   hh(rr) = plot(days_v(1):0.001:days_v(end), y_p', 'Color', colors(rr, :));
            else
                   plot(days_v, data_p', 'Color', colors(rr, :))
            end
            hold on
          
            label{rr} = ['Region ', num2str(regions(rr))];
        end
        title(type) 
        ylim([0.9 2.5])

    end
    
%     subplot(1, numel(types)+1, numel(types)+1)
    legend(hh, label, 'location', 'eastoutside')
    legend boxoff
    i = i+ n_v;
    
    suptitle(['Subject ' num2str(ss)])
    
    if rep
        fig = Figure(gcf);
        fig.Width = '8in';
        fig.Height = '5in';
        fig.Scaling = 'custom';
        add(rpt, fig);
    end
    
end

if rep
    close(rpt);
end