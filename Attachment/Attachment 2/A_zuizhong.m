clear all;
close all;

% 创建并打开一个用于保存结果的CSV文件，CSV文件可以方便地被电子表格软件打开并以表格形式查看
fileID = fopen('result.csv', 'w');

% 写入CSV文件的表头
fprintf(fileID, '文件名,是否可能偏色（均值差异）,是否可能偏色（方差比值）,是否偏色可能性较高（综合判断）,是否可能低亮度（均值判断）,是否可能低亮度（直方图判断）,偏色均值差异得分,偏色方差比值得分,偏色综合得分,弱光均值得分,弱光直方图得分,是否模糊\n');

% 设置文件夹路径，将这里替换为你实际存放图片的文件夹路径
folder_path = 'E:\附件一\附件'; 

% 获取文件夹中的所有.jpg和.png格式的图片文件列表
file_list = dir(fullfile(folder_path, '*.jpg'));
file_list = [file_list; dir(fullfile(folder_path, '*.png'))];

% 设定均值差值阈值和方差比值阈值
meanDiffThreshold = 45;
varRatioThreshold = 2;

% 设定判断弱光的阈值
lowLightThreshold = 80; 
histogramThreshold = 0.7; 

% 设定模糊判断的阈值
blurThreshold = 0.007; 

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
    meanDiffResult = '否';
    meanDiffScore = 0;
    if abs(redMean - greenMean) > meanDiffThreshold || abs(redMean - blueMean) > meanDiffThreshold || abs(greenMean - blueMean) > meanDiffThreshold
        meanDiffResult = '是';
        meanDiffScore = sum([abs(redMean - greenMean) / meanDiffThreshold, abs(redMean - blueMean) / meanDiffThreshold, abs(greenMean - blueMean) / meanDiffThreshold]);
    end

    % 判断方差比值是否超过阈值
    varRatioResult = '否';
    varRatioScore = 0;
    if redVar / greenVar > varRatioThreshold || redVar / blueVar > varRatioThreshold || greenVar / blueVar > varRatioThreshold
        varRatioResult = '是';
        varRatioScore = sum([redVar / greenVar / varRatioThreshold, redVar / blueVar / varRatioThreshold, greenVar / blueVar / varRatioThreshold]);
    end

    % 找到直方图峰值对应的颜色强度值
    [~, redPeakBin] = max(histRed);
    [~, greenPeakBin] = max(histGreen);
    [~, bluePeakBin] = max(histBlue);

    % 结合直方图峰值和统计差异综合判断偏色
    comprehensiveResult = '否';
    comprehensiveScore = 0;
    if (abs(redMean - greenMean) > meanDiffThreshold || abs(redMean - blueMean) > meanDiffThreshold || abs(greenMean - blueMean) > meanDiffThreshold) ||...
            (redVar / greenVar > varRatioThreshold || redVar / blueVar > varRatioThreshold || greenVar / blueVar > varRatioThreshold) ||...
            (abs(redPeakBin - greenPeakBin) > 30 || abs(redPeakBin - bluePeakBin) > 30 || abs(greenPeakBin - bluePeakBin) > 30)
        comprehensiveResult = '是';
        comprehensiveScore = meanDiffScore + varRatioScore + sum([abs(redPeakBin - greenPeakBin) / 30, abs(redPeakBin - bluePeakBin) / 30, abs(greenPeakBin - bluePeakBin) / 30]);
    end

    % 判断是否为弱光图像（针对海底图像相关处理部分）
    % 将彩色图像转换为灰度图像
    grayImage = rgb2gray(image);

    % 计算灰度图像像素值的均值
    meanGrayValue = mean(grayImage(:));

    % 计算灰度图像的直方图
    grayHist = imhist(grayImage);

    % 判断图像是否为弱光图像（基于均值）
    meanLightResult = '否';
    meanLightScore = 0;
    if meanGrayValue < lowLightThreshold
        meanLightResult = '是';
        meanLightScore = (lowLightThreshold - meanGrayValue) / lowLightThreshold;
    end

    % 计算灰度直方图中低亮度区域（这里假设0 - 50为低亮度区域）的像素比例
    lowLightPixelsRatio = sum(grayHist(1:50))/numel(grayImage);

    % 判断图像是否为弱光图像（基于直方图）
    histogramLightResult = '否';
    histogramLightScore = 0;
    if lowLightPixelsRatio > histogramThreshold
        histogramLightResult = '是';
        histogramLightScore = (lowLightPixelsRatio - histogramThreshold) / histogramThreshold;
    end

    % 判断图像是否模糊
    if size(image, 3) == 3
        img_gray = rgb2gray(image);
    else
        img_gray = image;
    end

    % 对图像进行傅里叶变换
    F = fft2(double(img_gray));
    F_shifted = fftshift(F);
    magnitude_spectrum = log(1 + abs(F_shifted));

    % 计算高频分量占比
    [rows, cols] = size(img_gray);
    center_row = floor(rows / 2);
    center_col = floor(cols / 2);
    radius = min(center_row, center_col) / 4;
    [X, Y] = meshgrid(1:cols, 1:rows);
    distance = sqrt((X - center_col).^2 + (Y - center_row).^2);

    % 高频区域（距离中心大于radius的部分）
    high_freq = distance > radius;
    high_freq_energy = sum(sum(abs(F_shifted(high_freq)).^2));
    total_energy = sum(sum(abs(F_shifted).^2));

    % 高频能量占比
    high_freq_ratio = high_freq_energy / total_energy;

    % 判断模糊程度
    blurResult = '否';
    if high_freq_ratio < blurThreshold
        blurResult = '是';
    end

    % 将当前图片的判断结果写入CSV文件
    fprintf(fileID, '%s,%s,%s,%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%s\n', file_name, meanDiffResult, varRatioResult, comprehensiveResult, meanLightResult, histogramLightResult, meanDiffScore, varRatioScore, comprehensiveScore, meanLightScore, histogramLightScore, blurResult);
end

% 关闭文件
fclose(fileID);