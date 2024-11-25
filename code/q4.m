clear;
clear all;
close all;
fileID = fopen('q4.csv', 'w');
fprintf(fileID, 'PSNR,UCIQE,UIQM\n');
% figure;
% size = [800,1200];
hold on;
% num = 3;

for i = 1:12
    fileName = sprintf('test_%03d.png', i);
    A = imread(fileName);
    A_1 = uwredcomp(A);
    % A_1 = imresize(A_1,size);
    % subplot(2,num,i)
    % imshow(A)
    % subplot(2,num,i+num)
    % imshow(A_1)
    img_original = A;
    img_processed = A_1;
    % figure;
    % imshow(img_original);
    % title('original image');
    %
    % figure;
    % imshow(img_processed);
    % title('processed image');
    img_original = im2double(img_original);
    img_processed = im2double(img_processed);
    uiqm_original = computeUIQM(img_original);
    uiqm_processed = computeUIQM(img_processed);
    psnr_value = computePSNR(img_original, img_processed);
    uciqe_original = computeUCIQE(img_original);
    uciqe_processed = computeUCIQE(img_processed);
    %fprintf('Quality assessment results\n');
    % fprintf('PSNR: %.2f dB\n\n', psnr_value);
    %fprintf('UCIQE of original image: %.4f\n', uciqe_original);
    % fprintf('UCIQE for processed images: %.4f\n', uciqe_processed);
    %fprintf('UIQM of original image: %.4f\n', uiqm_original);
    %fprintf('UIQM for processed images: %.4f\n\n', uiqm_processed);
    fprintf(fileID, '%.4f,%.4f,%.4f\n',psnr_value,uciqe_processed, uiqm_processed);
end
fclose(fileID);
% saveas(gcf,'q4.jpg')


function RGBnew = uwredcomp(RGB,varargin)
smoothing = 0;
kclip = 0.01;
gammaadj = 1;
alpha = 0.4;
if nargin>1
    k = 1;
    while k <= numel(varargin)
        thisarg = lower(varargin{k});
        switch thisarg
            case 'smoothing'
                smoothing = varargin{k+1};
                k = k+2;
            case 'kthresh'
                kclip = varargin{k+1};
                k = k+2;
            case 'gamma'
                gammaadj = varargin{k+1};
                k = k+2;
            case 'alpha'
                alpha = varargin{k+1};
                k = k+2;
            otherwise
                error('UWREDCOMP: unknown option %s',thisarg)
        end
    end
end
if smoothing ~= 0 && (ifversion('<','R2014a') || ~hasipt())
    error('UWREDCOMP: smoothing option requires imguidedfilter() from IPT (R2014a or newer)');
end
if size(RGB,3)~=3
    error('UWREDCOMP: expected INPICT to be RGB')
end
[RGB inclass] = imcast(RGB,'double');
lam = [620 540 450];
blam = (-0.00113*lam + 1.62517);
Blam = ctflop(quantile(RGB,1-0.001,[1 2]));
cgcr = (blam(2)*Blam(1))/(blam(1)*Blam(2));
cbcr = (blam(3)*Blam(1))/(blam(1)*Blam(3));
wrgb = [1 cgcr cbcr]/(1 + cgcr + cbcr);
Rnew = imappmat(RGB,wrgb);
if smoothing > 0
    Rnew = imguidedfilter(Rnew,RGB(:,:,2),'degree',smoothing);
end
RGBnew = RGB;
RGBnew(:,:,1) = Rnew;
inlim = stretchlimFB(RGBnew,kclip);
RGBcont = imadjustFB(RGBnew,inlim,[0 1],gammaadj);
RGBahq = zeros(size(RGBnew));
for c = 1:3
    RGBahq(:,:,c) = adapthisteqFB(RGBnew(:,:,c),'distribution','rayleigh');
end
RGBnew = RGBcont*alpha + RGBahq*(1-alpha);
RGBnew = imcast(RGBnew,inclass);
end

function uiqm = computeUIQM(img)
if size(img,3) ~= 3
    error('The input image must be an RGB image');
end
uicm = computeUICM(img);
uism = computeUISM(img);
uiconm = computeUIConM(img);
uiqm = 0.0282 * uicm + 0.2953 * uism + 3.5753 * uiconm;
end

function uicm = computeUICM(img)
R = img(:,:,1);
G = img(:,:,2);
B = img(:,:,3);
rg = R - G;
yb = 0.5 * (R + G) - B;
std_rg = std(rg(:));
std_yb = std(yb(:));
mean_rg = mean(rg(:));
mean_yb = mean(yb(:));
uicm = sqrt(std_rg^2 + std_yb^2) + 0.3 * sqrt(mean_rg^2 + mean_yb^2);
end

function uism = computeUISM(img)
sobel_x = [-1 0 1; -2 0 2; -1 0 1];
sobel_y = sobel_x';
gray = rgb2gray(img);
Ix = conv2(double(gray), sobel_x, 'same');
Iy = conv2(double(gray), sobel_y, 'same');
gradient_magnitude = sqrt(Ix.^2 + Iy.^2);
uism = mean(gradient_magnitude(:));
end

function uiconm = computeUIConM(img)
gray = rgb2gray(img);
uiconm = std(double(gray(:)));
end

function psnr_val = computePSNR(img1, img2)
if ~isequal(size(img1), size(img2))
    error('Both images must be the same size to calculate the PSNR.');
end
mse = mean( (img1(:) - img2(:)).^2 );
if mse == 0
    psnr_val = Inf;
else
    max_I = 1;
    psnr_val = 10 * log10( (max_I^2) / mse );
end
end

function uciqe = computeUCIQE(img)
if size(img,3) ~= 3
    error('The input image must be an RGB image');
end
lab = rgb2lab(img);
L = lab(:,:,1);
a = lab(:,:,2);
b = lab(:,:,3);
saturation = std([a(:), b(:)], 0, 1);
S = mean(saturation);
C = std(L(:));
sobel_x = [-1 0 1; -2 0 2; -1 0 1];
sobel_y = sobel_x';
Ix = conv2(double(L), sobel_x, 'same');
Iy = conv2(double(L), sobel_y, 'same');
gradient_magnitude = sqrt(Ix.^2 + Iy.^2);
H = mean(gradient_magnitude(:));
uciqe = 0.468 * S + 0.274 * C + 0.257 * H;
end
%imshow(dis)%,title('processed image')
