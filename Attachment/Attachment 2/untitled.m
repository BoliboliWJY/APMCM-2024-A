clear;

t = imread('test_001.png');  
% t = imread('image_080.png');  
I=t(:,:,1);
mean_red = mean(mean(I)) / 255;
Green = t(:,:,2);
mean_green = mean(mean(Green)) / 255;
Blue = t(:,:,3);
mean_blue = mean(mean(Blue)) / 255;
[height,width] = size(I);  

subplot(121);

imshow(t),title('原始图像')%显示原始图像   

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
figure

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

imshow(dis)%,title('处理之后的图像')%显示均衡化后的图像

