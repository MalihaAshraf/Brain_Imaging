%% Whole brain

types = {'AD', 'FA', 'MD', 'RD'};

n_visits  =  [16 14 14 13 12 13 13 14 13 13 13 11 12 11 10];
n_sub = length(n_visits);
load days.mat
days(2, :) = -2;
regions = [ 4 5 6 23 24 25 26 27 28 29 42 43 44 45];
colors = lines(n_sub);
rep = true;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Raw DATA

norm = false;
fit_ = false;

if rep
    import mlreportgen.report.*
    import mlreportgen.dom.*
    rpt = Report('Whole_brain','pdf');
    chap = Chapter('Raw data');
    add(rpt,chap);
    rpt.Document.CurrentPageLayout.PageMargins.Left = '0.5in';
    rpt.Document.CurrentPageLayout.PageMargins.Right = '0.5in';
end

figure('units','normalized','outerposition',[0 0 1 1]);

for tt = 1:numel(types)
    clf
    
    type = types{tt};
    for rr = 1:numel(regions)

        clear ds
        fileID = fopen([type '/' type '/ske', num2str(regions(rr)) '_' type '_values.txt'], 'r');
        formatSpec = '%f %f %f %f';
        sz = [4 Inf];
        ds = fscanf(fileID, formatSpec, sz);
        ds = ds';
        fclose(fileID);
        i = 1;
        subplot(3, 5, rr)
        axis square
        for ss = 1:n_sub
               n_v = n_visits(ss);
               data = ds(i+1:i+n_v-1, 1);
               data_sd = ds(i+1:i+n_v-1, 2);
               norm_param = data(2);
               
               if norm
                  data = (data-norm_param)./norm_param; 
               end
               
               days_v = days(:, ss);
               days_v = days_v(~isnan(days_v));
               days_v = days_v(2:end, :);
               
               scatter(days_v, data', 5, colors(ss, :), 'filled')
               hold on
               
               if fit_
                   ft = fittype('a/(1+exp(-b*x))');
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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Normalized DATA
if 1

norm = true;
fit_ = false;

if rep
    chap = Chapter('Normalized data');
    add(chap, Equation('$data = \frac{data - data_{ visit -1}}{data_{ visit -1}}$'));
    add(rpt,chap);
end

for tt = 1:numel(types)
    clf
    
    type = types{tt};
    for rr = 1:numel(regions)

        clear ds
        fileID = fopen([type '/' type '/ske', num2str(regions(rr)) '_' type '_values.txt'], 'r');
        formatSpec = '%f %f %f %f';
        sz = [4 Inf];
        ds = fscanf(fileID, formatSpec, sz);
        ds = ds';
        fclose(fileID);
        i = 1;
        subplot(3, 5, rr)
        axis square
        
        for ss = 1:n_sub
               n_v = n_visits(ss);
               data = ds(i+1:i+n_v-1, 1);
               data_sd = ds(i+1:i+n_v-1, 2);
               norm_param = data(2);
               
               if norm
                  data = (data-norm_param)./norm_param; 
               end
               
               days_v = days(:, ss);
               days_v = days_v(~isnan(days_v));
               days_v = days_v(2:end, :);
               
               scatter(days_v, data', 5, colors(ss, :), 'filled')
               hold on
               
               if fit_
                   ft = fittype('a/(1+exp(-b*x))');
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
    end
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Fitted DATA

if 1
norm = true;
fit_ = true;

if rep
    chap = Chapter('Fitted data');
    add(chap, Equation('$y = \frac{a}{1+e^{-bx}}$'));
    add(rpt,chap);
end

for tt = 1:numel(types)
    clf
    
    type = types{tt};
    for rr = 1:numel(regions)

        clear ds
        fileID = fopen([type '/' type '/ske', num2str(regions(rr)) '_' type '_values.txt'], 'r');
        formatSpec = '%f %f %f %f';
        sz = [4 Inf];
        ds = fscanf(fileID, formatSpec, sz);
        ds = ds';
        fclose(fileID);
        i = 1;
        subplot(3, 5, rr)
        axis square
        
        for ss = 1:n_sub
               n_v = n_visits(ss);
               data = ds(i+1:i+n_v-1, 1);
               data_sd = ds(i+1:i+n_v-1, 2);
               norm_param = data(2);
               
               if norm
                  data = (data-norm_param)./norm_param; 
               end
               
               days_v = days(:, ss);
               days_v = days_v(~isnan(days_v));
               days_v = days_v(2:end, :);
               
               scatter(days_v, data', 5, colors(ss, :), 'filled')
               hold on
               
               if fit_
                   ft = fittype('a/(1+exp(-b*x))');
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
    end
end

end

if rep
    close(rpt);
end

% 
%% Corpus Collasum


