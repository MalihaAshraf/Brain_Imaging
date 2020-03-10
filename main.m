%% Whole brain

types = {'AD', 'FA', 'MD', 'RD'};
n_visits  =  [16 14 14 13 12 13 13 14 13 13 13 11 12 11 10];
n_sub = length(n_visits);
load days.mat
regions = [ 4 5 6 23 24 25 26 27 28 29 42 43 44 45];
colors = lines(n_sub);
norm = true;
fit_ = true;

for tt = 1%:numel(types)
    figure('units','normalized','outerposition',[0 0 1 1])
    
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

        for ss = 1:n_sub
               n_v = n_visits(ss);
               data = ds(i:i+n_v-1, 1);
               data_sd = ds(i:i+n_v-1, 2);
               norm_param = data(2);
               
               if norm
                  data = (data-norm_param)./norm_param; 
               end
               
               days_v = days(:, ss);
               days_v = days_v(~isnan(days_v));

               scatter(days_v, data', 5, colors(ss, :), 'filled')
               hold on
               
               if fit_
                   f = fit((days_v), data, 'exp1', 'Weights', data_sd);
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
    
end



%% Corpus Collasum


