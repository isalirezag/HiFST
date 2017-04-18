% Here we use the functions which provided in tookbox at:
% https://github.com/yaoliUoA/evalsaliency

clear all
close all;
clc;

method = '';% e.g. 'HiFSTBlurMap';
dataset = 'BlurDetection';

GThpath = ['.\GT\*.png'];
savepath = ['.\output\'];
Result_path = ['..\Resutls\' '*.png'];


if ~exist(savepath,'dir')
    mkdir(savepath);
end
dir_im = dir(Result_path);
assert(~isempty(dir_im),'No  map found, please check the path!');
dir_tr= dir(GThpath);
assert(~isempty(dir_tr),'No ground-truth image found, please check the path!');
% assert(length(dir_im)==length(dir_tr),'The number of maps and ground-truth images are not equal!')
imNum = length(dir_tr);
precision = zeros(256,1);
recall = zeros(256,1);

%% compute PR curve
for i = 1:1000
    GTImgName = dir_tr(i).name;
    Map_im=[];
    GT_Img =  (imread([GThpath(1:end-5),GTImgName]));
%     figure,imshow(GT_Img,[])
    if length(size(GT_Img))==3
        GT_Img=rgb2gray(GT_Img);
    end
    %     Making sure the GT is binary
    GT_Img(GT_Img>0.5)=1;
    GT_Img(GT_Img~=1)=0;
    GT_Img= (GT_Img);
    if max(max(GT_Img))==255
        GT_Img = (GT_Img./255);
    end
    GT_Img=double(GT_Img);
    
    Map_im = (imread([Result_path(1:end-5), GTImgName(1:end-4),'_',method, '.png']));
    
    %     Making sure the Map is 2D (grayscale)
    if length(size(Map_im))==3
        Map_im=rgb2gray(Map_im);
    end
    
    Map_im = double(255-double(mat2gray( Map_im).*255));
    Map_im=imresize(Map_im,[size(GT_Img)]);
    
%     figure,imshow(Map_im,[])
%     return
    for threshold = 0:1:255
        index1 = (Map_im>=threshold);
        truePositive = length(find(index1 & GT_Img));
        groundTruth = length(find(GT_Img));
        detected = length(find(index1));
        if truePositive~=0
            precision(threshold+1) = precision(threshold+1)+truePositive/detected;
            recall(threshold+1) = recall(threshold+1)+truePositive/groundTruth;
        end
    end
    display(num2str(i));
end

precision = precision./imNum;
recall = recall./imNum;
pr = [precision'; recall'];
fid = fopen([savepath dataset,  method, '_PRCurve.txt'],'at');
fprintf(fid,'%f %f\n',pr);
fclose(fid);
disp('Done!');
