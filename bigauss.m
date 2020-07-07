function y = bigauss(a, b, c1, c2, x)
%UNTITLED6 Summary of this function goes here
%   Detailed explanation goes here
 
y1 = a.*exp(-((x-b)./c1).^2);
y2 = a.*exp(-((x-b)./c2).^2);
y = y1;
y(x>b) = y2(x>b);


end

