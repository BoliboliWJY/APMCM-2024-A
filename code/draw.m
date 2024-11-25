clc;
clear all;
close all;

f=imread('test_001.png');
subplot(211);
imshow(f);title('Image');
fcal=double(f); 
[m,n,h]=size(f);
Y=zeros(h,256);

for k=1:h
    for i=1:m
        for j=1:n
            Y(k,fcal(i,j,k)+1)=Y(k,fcal(i,j,k)+1)+1; %哪一灰度级出现一次其相应像素点数+1。灰度级范围是0-255，但矩阵是1-256，列数要额外+1
        end
    end
end

X=0:1:255; 
subplot(212);
histogram=bar(X,Y); 
axis([0 255,-inf inf]) 
% xlabel('gray level');ylabel('number of pixels');

if h==3 
    title('Color Histogram');
    set(histogram(1),'FaceColor',[1 0.1882 0.1882]); 
    set(histogram(2),'FaceColor',[0.5 1 0]);
    set(histogram(3),'FaceColor',[0 0.5 1]);
    FreqNum=zeros(size(f,3),256);

    for i=1:size(f,3)
        for j=0:255
            FreqNum(i,j+1)=sum(sum(f(:,:,i)==j));
        end
    end
    hold on
    plot(X,Y(1,:),'Color',[1 0.1882 0.1882]);  %加上边界轮廓
    plot(X,Y(2,:),'Color',[0.5 1 0]);
    plot(X,Y(3,:),'Color',[0 0.5 1]);
    hold off
else
    title('grey histogram');
end

