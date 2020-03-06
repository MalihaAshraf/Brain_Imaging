%%

n_visits  =  [16 14 14 13 12 13 13 14 13 13 13 11 12 11 10];
n_sub = length(n_visits);
regions = [ 4 5 6 23 24 25 26 27 28 29 42 43 44 45];

%%
clear ds
fileID = fopen('FA/FA/ske4_FA_values.txt', 'r');
formatSpec = '%f %f %f %f';
size = [4 Inf];
ds = fscanf(fileID, formatSpec, size);
ds = ds';
fclose(fileID);

%%
ds = dataset('file', 'FA/FA/ske4_FA_values.txt', 'Delimiter', ' ');