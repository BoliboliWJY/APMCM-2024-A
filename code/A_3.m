clear;
clear all;
close all;
fileID = fopen('iddddd.csv', 'w');
fprintf(fileID, 'PSNR,UCIQE,UIQM\n');
for ii=1:272
% for ij=273:400
fileName = sprintf('image_%03d.png',ii);
% fileName = sprintf('image_%03d.jpg',ij);
    t = imread(fileName);   
% t = imread('image_080.png');  
I=t(:,:,1);
mean_red = mean(mean(I)) / 255;
Green = t(:,:,2);
mean_green = mean(mean(Green)) / 255;
Blue = t(:,:,3);
mean_blue = mean(mean(Blue)) / 255;
[height,width] = size(I);  

%subplot(121);

%imshow(t),title('原始图像')%显示原始图像   

%对R通道进行均衡化处理，均衡化可以写一个统一的函数，直接调用

%进行像素灰度统计;  

s = zeros(1,256);%统计各灰度数目，共256个灰度级  
sg = zeros(1,256);
%绘制直方图

gp=zeros(1,256);

for k=0:255

    gp(k+1)=length(find(I==k))/(height*width);

end

for i = 1:height  

    for j = 1: width  

        s(I(i,j) + 1) = s(I(i,j) + 1) + 1;%对应灰度值像素点数量增加一  
        sg(Green(i,j) + 1) = sg(Green(i,j) + 1) + 1;
    end  

end  

%计算灰度分布密度  

p = zeros(1,256);  
pg = zeros(1,256);
for i = 1:256  

    p(i) = s(i) / (height * width * 1.0);  
    pg(i) = s(i) / (height * width * 1.0);      
end  

%计算累计直方图分布  

c = zeros(1,256);  
cg = zeros(1,256);
c(1) = p(1);

for i = 2:256   

        c(i) = c(i - 1) + p(i);  
        cg(i) = c(i - 1) + p(i);
end  

%累计分布取整,将其数值归一化为1~256 

c = uint8(255 .* c + 0.5);  
cg = uint8(255 .* c + 0.5);
%对图像进行均衡化
clear Ir Ig Ib dis;
for i = 1:height  

    for j = 1: width  

        Ir(i,j) = c(I(i,j)+1);  
        Ig(i,j) = cg(I(i,j)+1);
    end  

end  
%%
    Ir = Ir * (1 - mean_red) * 0.5 ; 
% Ir = Ir * 0.5; 

% Ig = Ig * 1;
dis(:,:,1)=Ir;
% dis(:,:,2)=Ig;

%subplot(122)  

%imshow(Ir)%显示均衡化后的图像

%对G通道进行均衡化处理，均衡化可以写一个统一的函数，直接调用

I=t(:,:,2);

[height,width] = size(I);  



%下面使用直方图均衡化进行处理

%进行像素灰度统计;  

s = zeros(1,256);%统计各灰度数目，共256个灰度级  

%绘制直方图

gp=zeros(1,256);

for k=0:255

    gp(k+1)=length(find(I==k))/(height*width);

end

for i = 1:height  

    for j = 1: width  

        s(I(i,j) + 1) = s(I(i,j) + 1) + 1;%对应灰度值像素点数量增加一  

    end  

end  

%计算灰度分布密度  

p = zeros(1,256);  

for i = 1:256  

    p(i) = s(i) / (height * width * 1.0);  

end  

%计算累计直方图分布  

c = zeros(1,256);  

c(1) = p(1);

for i = 2:256   

        c(i) = c(i - 1) + p(i);  

end  

%累计分布取整,将其数值归一化为1~256 

c = uint8(255 .* c + 0.5);  

%对图像进行均衡化

for i = 1:height  

    for j = 1: width  

        Ig(i,j) = c(I(i,j)+1);  

    end  

end  

%subplot(122)  

%imshow(Ig)%显示均衡化后的图像
%%
dis(:,:,2)=Ig * (1 - mean_green) * 0.5;



%对B通道进行均衡化处理，均衡化可以写一个统一的函数，直接调用

I=t(:,:,3);

[height,width] = size(I);  



%下面使用直方图均衡化进行处理

%进行像素灰度统计;  

s = zeros(1,256);%统计各灰度数目，共256个灰度级  

%绘制直方图

gp=zeros(1,256);

for k=0:255

    gp(k+1)=length(find(I==k))/(height*width);

end

for i = 1:height  

    for j = 1: width  

        s(I(i,j) + 1) = s(I(i,j) + 1) + 1;%对应灰度值像素点数量增加一  

    end  

end  

%计算灰度分布密度  

p = zeros(1,256);  

for i = 1:256  

    p(i) = s(i) / (height * width * 1.0);  

end  

%计算累计直方图分布  

c = zeros(1,256);  

c(1) = p(1);

for i = 2:256   

        c(i) = c(i - 1) + p(i);  

end  

%累计分布取整,将其数值归一化为1~256 

c = uint8(255 .* c + 0.5);  

%对图像进行均衡化
for i = 1:height  

    for j = 1: width  

        Ib(i,j) = c(I(i,j)+1);  

    end  

end  
%%
dis(:,:,3)=Ib * (1 - mean_blue) * 0.5; 

%subplot(122)  

%imshow(Ib)%显示均衡化后的图像 

% subplot(122);
%figure

[height, width, ch] = size(dis);
if ch == 3
    % 将RGB图像转换为双精度类型，方便后续计算
    image_double = im2double(dis);
    r_mean = mean(image_double(:,:,1));
    g_mean = mean(image_double(:,:,2));
    b_mean = mean(image_double(:,:,3));
    mean_grey = (r_mean + g_mean + b_mean)/3;
    image_double(:,:,1) = image_double(:,:,1) * (mean_grey/r_mean);
    image_double(:,:,2) = image_double(:,:,2) * (mean_grey/g_mean);
    image_double(:,:,3) = image_double(:,:,3) * (mean_grey/b_mean);
    dis = uint8(image_double*600);
end

% 读取原始图像和增强后图像
img_original = t;
img_processed = dis;

% 显示图像（可选）
% figure;
% imshow(img_original);
% title('原始图像');
% 
% figure;
% imshow(img_processed);
% title('增强后图像');

% 将图像转换为双精度并归一化到 [0,1]
img_original = im2double(img_original);
img_processed = im2double(img_processed);

% 计算 UIQM
uiqm_original = computeUIQM(img_original);
uiqm_processed = computeUIQM(img_processed);

% 计算 PSNR
psnr_value = computePSNR(img_original, img_processed);

% 计算 UCIQE
uciqe_original = computeUCIQE(img_original);
uciqe_processed = computeUCIQE(img_processed);

% 显示结果
%fprintf('--- 质量评估结果 ---\n');
% fprintf('原始与增强后图像的PSNR: %.2f dB\n\n', psnr_value);
%fprintf('原始图像的UCIQE: %.4f\n', uciqe_original);
% fprintf('增强后图像的UCIQE: %.4f\n', uciqe_processed);
%fprintf('原始图像的UIQM: %.4f\n', uiqm_original);
fprintf('增强后图像的UIQM: %.4f\n\n', uiqm_processed);
fprintf(fileID, '%.4f,%.4f,%.4f\n',psnr_value,uciqe_processed, uiqm_processed);
end
fclose(fileID);

function uiqm = computeUIQM(img)
    % 计算 UIQM（Underwater Image Quality Measure）
    % 确保图像是 RGB 图像
    if size(img,3) ~= 3
        error('输入图像必须是 RGB 图像');
    end
    
    % 计算颜色度 (UICM)
    uicm = computeUICM(img);
    
    % 计算清晰度 (UISM)
    uism = computeUISM(img);
    
    % 计算对比度 (UIConM)
    uiconm = computeUIConM(img);
    
    % UIQM 综合评分
    uiqm = 0.0282 * uicm + 0.2953 * uism + 3.5753 * uiconm;
end

function uicm = computeUICM(img)
    % 颜色度计算
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
    % 清晰度计算基于梯度
    % 使用 Sobel 算子计算梯度
    sobel_x = [-1 0 1; -2 0 2; -1 0 1];
    sobel_y = sobel_x';
    
    % 转为灰度图像
    gray = rgb2gray(img);
    
    Ix = conv2(double(gray), sobel_x, 'same');
    Iy = conv2(double(gray), sobel_y, 'same');
    
    gradient_magnitude = sqrt(Ix.^2 + Iy.^2);
    
    % 清晰度基于梯度均值
    uism = mean(gradient_magnitude(:));
end

function uiconm = computeUIConM(img)
    % 对比度计算
    % 使用标准差作为对比度度量
    gray = rgb2gray(img);
    uiconm = std(double(gray(:)));
end

function psnr_val = computePSNR(img1, img2)
    % 计算两幅图像的 PSNR（Peak Signal-to-Noise Ratio）
    % 确保两幅图像大小相同
    if ~isequal(size(img1), size(img2))
        error('两幅图像的大小必须相同才能计算 PSNR');
    end
    
    mse = mean( (img1(:) - img2(:)).^2 );
    
    if mse == 0
        psnr_val = Inf;
    else
        max_I = 1; % 因为图像已经归一化到 [0,1]
        psnr_val = 10 * log10( (max_I^2) / mse );
    end
end

function uciqe = computeUCIQE(img)
    % 计算 UCIQE（Underwater Color Image Quality Evaluation）
    % 无参考图像质量评价指标
    
    % 确保图像是 RGB 图像
    if size(img,3) ~= 3
        error('输入图像必须是 RGB 图像');
    end
    
    % 将 RGB 图像转换为 CIELAB 色彩空间
    lab = rgb2lab(img);
    
    % 提取 L, a, b 通道
    L = lab(:,:,1);
    a = lab(:,:,2);
    b = lab(:,:,3);
    
    % 计算 Saturation（色彩饱和度）: a 和 b 通道的标准差
    saturation = std([a(:), b(:)], 0, 1);
    S = mean(saturation);
    
    % 计算 Contrast（对比度）: L 通道的标准差
    C = std(L(:));
    
    % 计算 Sharpness（清晰度）: L 通道梯度幅值的均值
    % 使用 Sobel 算子计算梯度
    sobel_x = [-1 0 1; -2 0 2; -1 0 1];
    sobel_y = sobel_x';
    
    Ix = conv2(double(L), sobel_x, 'same');
    Iy = conv2(double(L), sobel_y, 'same');
    
    gradient_magnitude = sqrt(Ix.^2 + Iy.^2);
    H = mean(gradient_magnitude(:));
    
    % 计算 UCIQE
    uciqe = 0.468 * S + 0.274 * C + 0.257 * H;
end

%imshow(dis)%,title('处理之后的图像')%显示均衡化后的图像
