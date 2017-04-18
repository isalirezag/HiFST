% This code is written by: S. Alireza Golestaneh;
% Email: sgolest1@asu.edu;
% March 1 2017; for CVPR 2017; Paper: Spatially-Varying  Blur Detection Based
% on Multiscale Fused and Sorted Transform Coefficients of Gradient Magnitudes
% S. Alireza Golestaneh and Lina Karam

clc; clear; close all; path(pathdef)

ImgName='motion0138.jpg';
InputImg=imread(ImgName);
 
[InputImg FinalMap] = HiFST(InputImg,1,0);
figure,imshow(InputImg,[])
figure,imshow(FinalMap,[])

% imwrite(mat2gray(FinalMap),[ImgName(1:end-4) '_FinalMap.png'],'Mode','lossless' );

% Important Note
% For input ImgName='motion0138.jpg'; you must get the 
% mean2(FinalMap)=0.1809; otherwise, there is something
% different in your setting than the setting that we
% used for our experiment/analysis

