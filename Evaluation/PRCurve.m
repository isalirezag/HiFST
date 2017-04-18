% Here we use the functions which provided in tookbox at:
% https://github.com/yaoliUoA/evalsaliency

clc
clear all
close all
warning off
path(pathdef)


dataset = 'BlurDetection'; % name of the dataset


methods = {'HiFSTBlurMap','Chang'   ,'Yi','JNB', 'Shi','Tang','Zhuo','Su', 'Chakrabarti','Liu','Bae' }; % you can add more names of methods separated by comma
methods2 = {'HiFST (Proposed)' , 'Tang [32]','Yi [36]' ,'Shi [26]','Shi [27]','Tang [31]','Zhuo [43]','Su [28]','Chakrabarti [3]','Liu [20]','Bae [1]'}; % you can add more names of methods separated by comma


methods_colors = distinguishable_colors(length(methods));
readpath = '.\output\';

%% load PRCurve.txt and draw PR curves
figure
hold on
for m = 1:length(methods)
    prFileName = strcat(readpath,dataset, '_', methods{m}, '_PRCurve.txt');
    
    R = load(prFileName);
    precision = R(:, 1);
    recall = R(:, 2);
    plot(recall, precision,'color',methods_colors(m,:),'linewidth',2);
end
axis([0 1 0 0.92]);
hold off
grid on;

legend( methods2 , 'Location', 'SouthWest');
title('Precision-Recall','fontsize',14,'fontname','times')
xlabel('Recall','fontsize',14,'fontname','times');
ylabel('Precision','fontsize',14,'fontname','times');
