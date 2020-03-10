
%%

types = {'AD', 'FA', 'MD', 'RD'};
n_visits  =  [16 14 14 13 12 13 13 14 13 13 13 11 12 11 10];
n_sub = length(n_visits);
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
        
       
        for rr = 1:numel(regions)
            fileID = fopen([type '/' type '/ske', num2str(regions(rr)) '_' type '_values.txt'], 'r');
            formatSpec = '%f %f %f %f';
            size = [4 Inf];
            ds = fscanf(fileID, formatSpec, size);
            ds = ds';
            fclose(fileID);
            
            data = ds(i:i+n_v-1, 1);
            data_sd = ds(i:i+n_v-1, 2);
            norm_param = data(1);
            
            if norm
                data = (data-norm_param)./norm_param; 
                data = data+0.1*(rr-1);
            end
            
            scatter(1:n_v, data', 5, colors(rr, :), 'filled')
            hold on
            
            if fit_
                   f = fit((1:n_v)', data, 'exp1', 'Weights', data_sd);
                   y = feval(f, 1:0.1:n_v);
                   hh(rr) = plot(1:0.1:n_v, y', 'Color', colors(rr, :));
            else
                   plot(1:n_v, data', 'Color', colors(rr, :))
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
regions = [ 4 5 6 23 24 25 26 27 28 29 42 43 44 45];
colors = lines(numel(regions));
norm = true;
fit_ = true;


i = 1;
for ss = 1:5
    
    figure('units','normalized','outerposition',[0 0 1 1])

    n_v = n_visits(ss);
    
    for tt = 1:numel(types)
        clear ds
%         subplot(1, n_sub, ss+(tt-1)*n_sub)
        subplot(1, numel(types), tt)
         type = types{tt};
         
        for rr = 1:numel(regions)
            fileID = fopen([type '/' type '/ske', num2str(regions(rr)) '_' type '_values.txt'], 'r');
            formatSpec = '%f %f %f %f';
            size = [4 Inf];
            ds = fscanf(fileID, formatSpec, size);
            ds = ds';
            fclose(fileID);
            
            data = ds(i:i+n_v-1, 1);
            data_sd = ds(i:i+n_v-1, 2);
            norm_param = data(1);
            
            if norm
                data = (data-norm_param)./norm_param; 
                data = data+0.1*(rr-1);
            end
            
            scatter(1:n_v, data', 5, colors(rr, :), 'filled')
            hold on
            
            if fit_
                   f = fit((1:n_v)', data, 'exp1', 'Weights', data_sd);
                   y = feval(f, 1:0.1:n_v);
                   hh(rr) = plot(1:0.1:n_v, y', 'Color', colors(rr, :));
            else
                   plot(1:n_v, data', 'Color', colors(rr, :))
            end
            hold on
          
            label{rr} = ['Region ', num2str(regions(rr))];
        end
        title(type) 
        ylim([-0.2 1.6])

    end
    legend(hh, label)
    legend boxoff
    i = i+ n_v;
    
    suptitle(['Sub ' num2str(ss)])
end