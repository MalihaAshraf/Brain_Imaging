function [y_out, mean_y, sign_y, peak_y, min_y] = process_y_bigauss(y)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

mean_y = mean(y);

n = 3;
y_ = y(1:n);
mean_y_ = mean_y;
neg_mean = abs(mean(y_ (y_ < mean_y_)-mean_y_));
pos_mean = abs(mean(y_ (y_ > mean_y_)-mean_y_));

if neg_mean > pos_mean
    sign_y = -1;
    peak_y = neg_mean;
else
    sign_y = 1;
    peak_y = pos_mean;
end

y_out = y-mean_y;
y_out = y_out.*sign_y;

min_y = min(y_out);
y_out = y_out - min_y;


end

