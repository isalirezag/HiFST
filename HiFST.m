% This code is written by: S. Alireza Golestaneh;
% Email: sgolest1@asu.edu;
% March 1 2017; for CVPR 2017; Paper: Spatially-Varying  Blur Detection Based
% on Multiscale Fused and Sorted Transform Coefficients of Gradient Magnitudes
% S. Alireza Golestaneh and Lina Karam

function [InputImg FinalMap] =...
    HiFST(InputImg,Slid,Noise)
% You can put [InputImg FinalMap  T D Camera_Point Camera_Point_OnImg
% Score] for the out put of the HiFST function instead of [InputImg FinalMap]
% to have more outputs if you like

InputImg = imnoise(InputImg,'gaussian',0,Noise);

InputImgGaus=(imgaussfilt(InputImg,0.5));
if length(size(InputImgGaus))==3
    InputImgGray=rgb2gray(InputImgGaus);
else
    InputImgGray=InputImgGaus;
end
% Computing Gradient Magnitude
GMImg = imgradient(double(InputImgGray));

% Change Flag to 0, if interested in using power of 2 block size.
% Please note that the results in the paper are produced by Flag=1;
Flag=1;
if Flag==1
    M_1=(2^3)-1;    M_2=(2^4)-1;    M_3=(2^5)-1;    M_4=(2^6)-1;
else Flag=0
    M_1=(2^3);    M_2=(2^4);    M_3=(2^5);    M_4=(2^6);
end
SelectedNumLayers=M_1+M_2+M_3+M_4+1;

% Labling the coeffitients to select the high frequencies
OutIndex1=FreqBands(M_1);
OutIndex2=FreqBands(M_2);
OutIndex3=FreqBands(M_3);
OutIndex4=FreqBands(M_4);

GMImg_P=padarray(GMImg,[floor(M_4/2) floor(M_4/2)]);
Padded_Size=size(GMImg_P);

n=0;
for i= floor(max([M_1,M_2,M_3,M_4])/2)+1:Slid:Padded_Size(1)-floor(max([M_1,M_2,M_3,M_4])/2)
    m=0; n=n+1;
    %        Showing the process
    disp_progress; clc; disp_progress(i, Padded_Size(1));
    
    for j= floor(max([M_1,M_2,M_3,M_4])/2)+1:Slid:Padded_Size(2)-floor(max([M_1,M_2,M_3,M_4])/2)
        m=m+1;
        
        %         Block selection
        Patch1=(GMImg_P(i-floor(M_1/2):i+floor(M_1/2),j- floor(M_1/2):j+ floor(M_1/2)));
        Patch2=(GMImg_P(i-floor(M_2/2):i+floor(M_2/2),j- floor(M_2/2):j+ floor(M_2/2)));
        Patch3=(GMImg_P(i-floor(M_3/2):i+floor(M_3/2),j- floor(M_3/2):j+ floor(M_3/2)));
        Patch4=(GMImg_P(i-floor(M_4/2):i+floor(M_4/2),j- floor(M_4/2):j+ floor(M_4/2)));
        
        %        Computing DCTs
        DCT_Coef1=abs(DCTA(Patch1));
        H1= DCT_Coef1( OutIndex1==0);
        
        DCT_Coef2=abs(DCTA(Patch2));
        H2= DCT_Coef2( OutIndex2==0);
        
        DCT_Coef3=abs((DCTA(Patch3)));
        H3= DCT_Coef3( OutIndex3==0);
        
        DCT_Coef4=abs((DCTA(Patch4)));
        H4= DCT_Coef4( OutIndex4==0);
        
        %       Sorting
        H_Sorted=sort([ H1 ; H2 ;H3 ;H4 ]);
        L(n,m,:) = H_Sorted(1:SelectedNumLayers);
        
    end
end

% Normalizing each layer, we only consider SelectedNumLayers number of
% layers
L_hat=[];
for i2=1:SelectedNumLayers
    L_hat(:,:,i2)=   mat2gray(L(:,:,i2)) ;
end

% Maxpooling
T = []; T = max(L_hat ,[],3);

%  Final Map and Post-Processing
D = [];
D = (entropyfilt(mat2gray(T),true(7))).*mat2gray(T);
Score=0;
Score=[median(mat2gray(D(:)))];

FinalMap=[];
FinalMap = RF(  (D),15, 0.25 , 3,imfilter(imresize( mat2gray(double(InputImgGray)),size(D)), fspecial('gaussian', [3, 3], 1)));

% Camera focus point map estiomation
P_Map=mat2gray(imgaussfilt(D,10));
Thr=0.98;
P_Map(P_Map>=Thr)=1;P_Map(P_Map~=1)=0;

% Binary Map
Camera_Point=mat2gray(P_Map);

% Camera focus point map estiomation on the Original Image
Camera_Point_OnImg=[];
A4= (imresize(mat2gray(InputImg(:,:,1)),size(Camera_Point)).*(1-Camera_Point));  A4( Camera_Point ==1) =1;
Camera_Point_OnImg(:,:,1)=A4;
Camera_Point_OnImg(:,:,2)=imresize(mat2gray(InputImg(:,:,2)),size(Camera_Point)).*(1-Camera_Point);
Camera_Point_OnImg(:,:,3)=imresize(mat2gray(InputImg(:,:,3)),size(Camera_Point)).*(1-Camera_Point);
end



%% Functions
function OutIndex=FreqBands(MatrixSize)
NewMatrixIndx=zeros(MatrixSize);
for i=1:MatrixSize
    NewMatrixIndx(1:((MatrixSize-1)/2)-i+2,i)=1;
end
for i=1:MatrixSize
    if (MatrixSize-((MatrixSize-1)/2)-i+01+0)<=0;
        NewMatrixIndx(1:MatrixSize-i+0,i)=2 ;
    else
        NewMatrixIndx(MatrixSize-((MatrixSize-1)/2)-i+1+0:MatrixSize-i+0,i)=2;
    end
end
NewMatrixIndx(1,1)=3;OutIndex=NewMatrixIndx;
end
function DCT_COEF=DCTA(Img_Blk1)
D=dctmtx(size(Img_Blk1,1));
DCT_COEF =  ( D*Img_Blk1*D');
end

%% The following function are obtain from the other peaople's code,
% which referenced in the paper.

function disp_progress(p, p_max)

persistent p_last;
if (nargin == 0)
    p_last = [];
    fprintf(1, '%s\n', '');
    return;
end

p_done = p / p_max * 100;
p_done = round(p_done / 10) * 10;

%[p_done p_last]

if (p_done == p_last)
    return;
end

if (~isempty(p_last))
    %  fprintf(1, '%d\n', p_last);
    fprintf(1, '%s\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b', '');
    %  return;
end
p_last = p_done;

switch (p_done)
    case 0
        fprintf(1, '%s', '[          ] 0%  ');
    case 10
        fprintf(1, '%s', '[|         ] 10% ');
    case 20
        fprintf(1, '%s', '[||        ] 20% ');
    case 30
        fprintf(1, '%s', '[|||       ] 30% ');
    case 40
        fprintf(1, '%s', '[||||      ] 40% ');
    case 50
        fprintf(1, '%s', '[|||||     ] 50% ');
    case 60
        fprintf(1, '%s', '[||||||    ] 60% ');
    case 70
        fprintf(1, '%s', '[|||||||   ] 70% ');
    case 80
        fprintf(1, '%s', '[||||||||  ] 80% ');
    case 90
        fprintf(1, '%s', '[||||||||| ] 90% ');
    case 100
        fprintf(1, '%s', '[||||||||||] 100% ');
end
drawnow;
end

%  RF  Domain transform recursive edge-preserving filter.
%
%  F = RF(img, sigma_s, sigma_r, num_iterations, joint_image)
%
%  Parameters:
%    img             Input image to be filtered.
%    sigma_s         Filter spatial standard deviation.
%    sigma_r         Filter range standard deviation.
%    num_iterations  Number of iterations to perform (default: 3).
%    joint_image     Optional image for joint filtering.
%
%
%
%  This is the reference implementation of the domain transform RF filter
%  described in the paper:
%
%    Domain Transform for Edge-Aware Image and Video Processing
%    Eduardo S. L. Gastal  and  Manuel M. Oliveira
%    ACM Transactions on Graphics. Volume 30 (2011), Number 4.
%    Proceedings of SIGGRAPH 2011, Article 69.
%
%  Please refer to the publication above if you use this software. For an
%  up-to-date version go to:
%
%             http://inf.ufrgs.br/~eslgastal/DomainTransform/
%
%
%  THIS SOFTWARE IS PROVIDED "AS IS" WITHOUT ANY EXPRESSED OR IMPLIED WARRANTIES
%  OF ANY KIND, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
%  FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL THE
%  AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
%  LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
%  OUT OF OR IN CONNECTION WITH THIS SOFTWARE OR THE USE OR OTHER DEALINGS IN
%  THIS SOFTWARE.
%
%  Version 1.0 - August 2011.

function F = RF(img, sigma_s, sigma_r, num_iterations, joint_image)

I = double(img);

if ~exist('num_iterations', 'var')
    num_iterations = 3;
end

if exist('joint_image', 'var') && ~isempty(joint_image)
    J = double(joint_image);
    
    if (size(I,1) ~= size(J,1)) || (size(I,2) ~= size(J,2))
        error('Input and joint images must have equal width and height.');
    end
else
    J = I;
end

[h w num_joint_channels] = size(J);

% Compute the domain transform (Equation 11 of our paper).
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Estimate horizontal and vertical partial derivatives using finite
% differences.
dIcdx = diff(J, 1, 2);
dIcdy = diff(J, 1, 1);

dIdx = zeros(h,w);
dIdy = zeros(h,w);

% Compute the l1-norm distance of neighbor pixels.
for c = 1:num_joint_channels
    dIdx(:,2:end) = dIdx(:,2:end) + abs( dIcdx(:,:,c) );
    dIdy(2:end,:) = dIdy(2:end,:) + abs( dIcdy(:,:,c) );
end

% Compute the derivatives of the horizontal and vertical domain transforms.
dHdx = (1 + sigma_s/sigma_r * dIdx);
dVdy = (1 + sigma_s/sigma_r * dIdy);

% We do not integrate the domain transforms since our recursive filter
% uses the derivatives directly.
%ct_H = cumsum(dHdx, 2);
%ct_V = cumsum(dVdy, 1);

% The vertical pass is performed using a transposed image.
dVdy = dVdy';

%% Perform the filtering.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N = num_iterations;
F = I;

sigma_H = sigma_s;

for i = 0:num_iterations - 1
    
    % Compute the sigma value for this iteration (Equation 14 of our paper).
    sigma_H_i = sigma_H * sqrt(3) * 2^(N - (i + 1)) / sqrt(4^N - 1);
    
    F = TransformedDomainRecursiveFilter_Horizontal(F, dHdx, sigma_H_i);
    F = image_transpose(F);
    
    F = TransformedDomainRecursiveFilter_Horizontal(F, dVdy, sigma_H_i);
    F = image_transpose(F);
    
end

F = cast(F, class(img));

end

% Recursive filter.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function F = TransformedDomainRecursiveFilter_Horizontal(I, D, sigma)

% Feedback coefficient (Appendix of our paper).
a = exp(-sqrt(2) / sigma);

F = I;
V = a.^D;

[h w num_channels] = size(I);

% Left -> Right filter.
for i = 2:w
    for c = 1:num_channels
        F(:,i,c) = F(:,i,c) + V(:,i) .* ( F(:,i - 1,c) - F(:,i,c) );
    end
end

% Right -> Left filter.
for i = w-1:-1:1
    for c = 1:num_channels
        F(:,i,c) = F(:,i,c) + V(:,i+1) .* ( F(:,i + 1,c) - F(:,i,c) );
    end
end

end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function T = image_transpose(I)

[h w num_channels] = size(I);

T = zeros([w h num_channels], class(I));

for c = 1:num_channels
    T(:,:,c) = I(:,:,c)';
end

end

%%
