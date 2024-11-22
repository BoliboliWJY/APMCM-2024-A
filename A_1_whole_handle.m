clear all;
close all;

% 设置文件夹路径，将这里替换为你实际存放图片的文件夹路径
folder_path = 'C:\Users\ASUS\Desktop\git\APMCM-2024-A\Attachment\Attachment 1'; 

% 获取文件夹中的所有.jpg和.png格式的图片文件列表
file_list = dir(fullfile(folder_path, '*.jpg'));
file_list = [file_list; dir(fullfile(folder_path, '*.png'))];

% 设定均值差值阈值和方差比值阈值
meanDiffThreshold = 20;
varRatioThreshold = 2;

% 设定判断弱光的阈值
lowLightThreshold = 50; 
histogramThreshold = 0.7; 

for i = 1:length(file_list)
    % 获取当前图片的完整文件名
    file_name = fullfile(folder_path, file_list(i).name);

    % 读取图像
    image = imread(file_name);

    % 数据类型转换（如果需要）
    if isa(image, 'uint16')
        image = im2uint8(image);
    end

    % 分离颜色通道
    redChannel = image(:, :, 1);
    greenChannel = image(:, :, 2);
    blueChannel = image(:, :, 3);

    % 计算颜色通道直方图
    histRed = imhist(redChannel);
    histGreen = imhist(greenChannel);
    histBlue = imhist(blueChannel);

    % 计算各通道均值
    redMean = mean(redChannel(:));
    greenMean = mean(greenChannel(:));
    blueMean = mean(blueChannel(:));

    % 计算各通道方差
    redVar = var(double(redChannel(:)), 1);
    greenVar = var(double(greenChannel(:)), 1);
    blueVar = var(double(blueChannel(:)), 1);

    % 判断均值差异是否超过阈值
    if abs(redMean - greenMean) > meanDiffThreshold || abs(redMean - blueMean) > meanDiffThreshold || abs(greenMean - blueMean) > meanDiffThreshold
        disp([file_name,'图像可能存在偏色（基于均值差异判断）']);
    end

    % 判断方差比值是否超过阈值
    if redVar / greenVar > varRatioThreshold || redVar / blueVar > varRatioThreshold || greenVar / blueVar > varRatioThreshold
        disp([file_name,'图像可能存在偏色（基于方差比值判断）']);
    end

    % 找到直方图峰值对应的颜色强度值
    [~, redPeakBin] = max(histRed);
    [~, greenPeakBin] = max(histGreen);
    [~, bluePeakBin] = max(histBlue);

    % 结合直方图峰值和统计差异综合判断偏色
    if (abs(redMean - greenMean) > meanDiffThreshold || abs(redMean - blueMean) > meanDiffThreshold || abs(greenMean - blueMean) > meanDiffThreshold) ||...
            (redVar / greenVar > varRatioThreshold || redVar / blueVar > varRatioThreshold || greenVar / blueVar > varRatioThreshold) ||...
            (abs(redPeakBin - greenPeakBin) > 30 || abs(redPeakBin - bluePeakBin) > 30 || abs(greenPeakBin - bluePeakBin) > 30)
        disp([file_name,'图像存在偏色可能性较高，综合判断依据：均值差异、方差比值及直方图峰值差异']);
    end
    
    % % 可视化直方图（可选，若要可视化需取消注释以下代码并调整figure相关设置以便正确显示）
    % figure;
    % subplot(3, 1, 1);
    % bar(histRed);
    % title('红色通道直方图');
    % subplot(3, 1, 2);
    % bar(histGreen);
    % title('绿色通道直方图');
    % subplot(3, 1, 3);
    % bar(histBlue);
    % title('蓝色通道直方图');

    % 判断是否为弱光图像（针对海底图像相关处理部分）
    % 将彩色图像转换为灰度图像
    grayImage = rgb2gray(image);

    % 计算灰度图像像素值的均值
    meanGrayValue = mean(grayImage(:));

    % 计算灰度图像的直方图
    grayHist = imhist(grayImage);

    % 判断图像是否为弱光图像（基于均值）
    if meanGrayValue < lowLightThreshold
        disp([file_name,'图像可能是弱光图像（基于均值判断）']);
    end

    % 计算灰度直方图中低亮度区域（这里假设0 - 50为低亮度区域）的像素比例
    lowLightPixelsRatio = sum(grayHist(1:50))/numel(grayImage);

    % 判断图像是否为弱光图像（基于直方图）
    if lowLightPixelsRatio > histogramThreshold
        disp([file_name,'图像可能是弱光图像（基于直方图判断）']);
    end
end