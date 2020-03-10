clear
figure
hold on 

for sub = 1:15 
clear
figure
load('S03MD_CC.mat')
plot(meanFADifValue,'DisplayName','meanFADifValue')
title ('S03MD_CC')

end  

clear
figure
load('S15MD_CC.mat')
plot(meanFADifValue,'DisplayName','meanFADifValue')
title ('S15MD_CC')