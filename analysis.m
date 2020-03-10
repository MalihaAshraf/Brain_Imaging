%% Dataset for PCA

types = {'AD', 'FA', 'MD', 'RD'};
n_visits  =  [16 14 14 13 12 13 13 14 13 13 13 11 12 11 10];
n_sub = length(n_visits);
regions = [ 4 5 6 23 24 25 26 27 28 29 42 43 44 45];
colors = lines(numel(regions));
norm = true;
fit_ = true;

dt = dataset();

i = 1;
k = 1;

for ss = 1:n_sub
    n_v = n_visits(ss);
    
    for tt = 1:numel(types)
        clear ds
        type = types{tt};
         
        for rr = 1:numel(regions)
            fileID = fopen([type '/' type '/ske', num2str(regions(rr)) '_' type '_values.txt'], 'r');
            formatSpec = '%f %f %f %f';
            sz = [4 Inf];
            ds = fscanf(fileID, formatSpec, sz);
            ds = ds';
            fclose(fileID);
            
            data = ds(i:i+n_v-1, 1);
            data_sd = ds(i:i+n_v-1, 2);
            norm_param = data(1);
            
            if norm
                data = (data-norm_param)./norm_param; 
            end
            
            if fit_
                f = fit((1:n_v)', data, 'exp1', 'Weights', data_sd);
                y = coeffvalues(f);
                dt.([type '_' num2str(regions(rr)), '_a'])(ss, 1) = y(1);
                dt.([type '_' num2str(regions(rr)), '_b'])(ss, 1) = y(2);
                D_mat(ss, tt, rr, 1) = y(1);
                D_mat(ss, tt, rr, 2) = y(2);
                var_ind(k, :) = [tt, rr, 1];
                var_ind(k+1, :) = [tt, rr, 2];
            end  
            k = k+2;
        end
    end
    i = i+ n_v;
end


%%

D = dataset2cell(dt);
vars = D(1, :);
var_ind = var_ind(1:length(vars), :);

D = cell2mat(D(2:end, :));
% [coeff,score,latent,tsquared,explained,mu] = pca(D);
[rho,pval] = corr(D);
sig_p = (pval < 0.05) & (abs(rho) > 0.85);
sig_p = tril(sig_p);

[i, j] = ind2sub(size(sig_p), find(sig_p));


%% Plot significant relation

figure,
fit_vars = {'a', 'b'};
for pp = 1:length(i)
    subplot(4, 6, pp)
    x = D_mat(:, var_ind(i(pp), 1),var_ind(i(pp), 2), var_ind(i(pp), 3));
    y = D_mat(:, var_ind(j(pp), 1),var_ind(j(pp), 2), var_ind(j(pp), 3));
    scatter(x, y)
    xlabel([types{var_ind(i(pp), 1)}, ' ', num2str(regions(var_ind(i(pp), 2))), ' ', fit_vars{ var_ind(i(pp), 3)}])
    ylabel([types{var_ind(j(pp), 1)}, ' ', num2str(regions(var_ind(j(pp), 2))), ' ', fit_vars{ var_ind(j(pp), 3)}])
end


