clear all;
close all;

fileID = fopen('result.csv', 'w');
fprintf(fileID, '文件名,是否可能偏色（均值差异）,是否可能偏色（方差比值）,是否偏色可能性较高（综合判断）,是否可能低亮度（均值判断）,是否可能低亮度（直方图判断）,偏色均值差异得分,偏色方差比值得分,偏色综合得分,弱光均值得分,弱光直方图得分,弱光最终得分,是否模糊,模糊程度分数,受损类别\n');
folder_path = 'C:\Users\22477\Desktop\Attachment 1'; 
file_list = dir(fullfile(folder_path, '*.jpg'));
file_list = [file_list; dir(fullfile(folder_path, '*.png'))];

meanDiffThreshold = 45;
varRatioThreshold = 2;

lowLightThreshold = 95; 
histogramThreshold = 0.7; 

blurThreshold = 0.003; 

max_d = 0;
for i = 1:length(file_list)
    file_name = fullfile(folder_path, file_list(i).name);
    image = imread(file_name);
    if isa(image, 'uint16')
        image = im2uint8(image);
    end

    redChannel = image(:, :, 1);
    greenChannel = image(:, :, 2);
    blueChannel = image(:, :, 3);

    histRed = imhist(redChannel);
    histGreen = imhist(greenChannel);
    histBlue = imhist(blueChannel);

    redMean = mean(redChannel(:));
    greenMean = mean(greenChannel(:));
    blueMean = mean(blueChannel(:));

    redVar = var(double(redChannel(:)), 1);
    greenVar = var(double(greenChannel(:)), 1);
    blueVar = var(double(blueChannel(:)), 1);

    meanDiffResult = '否';
    meanDiffScore = 0;
    if abs(redMean - greenMean) > meanDiffThreshold || abs(redMean - blueMean) > meanDiffThreshold || abs(greenMean - blueMean) > meanDiffThreshold
        meanDiffResult = '是';
        meanDiffScore = sum([abs(redMean - greenMean) / meanDiffThreshold, abs(redMean - blueMean) / meanDiffThreshold, abs(greenMean - blueMean) / meanDiffThreshold]);
    end

    varRatioResult = '否';
    varRatioScore = 0;
    if redVar / greenVar > varRatioThreshold || redVar / blueVar > varRatioThreshold || greenVar / blueVar > varRatioThreshold
        varRatioResult = '是';
        varRatioScore = sum([redVar / greenVar / varRatioThreshold, redVar / blueVar / varRatioThreshold, greenVar / blueVar / varRatioThreshold]);
    end

    [~, redPeakBin] = max(histRed);
    [~, greenPeakBin] = max(histGreen);
    [~, bluePeakBin] = max(histBlue);

    comprehensiveResult = '否';
    comprehensiveScore = 0;
    if (abs(redMean - greenMean) > meanDiffThreshold || abs(redMean - blueMean) > meanDiffThreshold || abs(greenMean - blueMean) > meanDiffThreshold) ||...
            (redVar / greenVar > varRatioThreshold || redVar / blueVar > varRatioThreshold || greenVar / blueVar > varRatioThreshold) ||...
            (abs(redPeakBin - greenPeakBin) > 30 || abs(redPeakBin - bluePeakBin) > 30 || abs(greenPeakBin - bluePeakBin) > 30)
        comprehensiveResult = '是';
        comprehensiveScore = meanDiffScore + varRatioScore + sum([abs(redPeakBin - greenPeakBin) / 30, abs(redPeakBin - bluePeakBin) / 30, abs(greenPeakBin - bluePeakBin) / 30]);
    end

    grayImage = rgb2gray(image);

    meanGrayValue = mean(grayImage(:));

    grayHist = imhist(grayImage);

    meanLightResult = '否';
    meanLightScore = 0;
    if meanGrayValue < lowLightThreshold
        meanLightResult = '是';
        meanLightScore = (lowLightThreshold - meanGrayValue) / lowLightThreshold;
    end

    lowLightPixelsRatio = sum(grayHist(1:50))/numel(grayImage);
    histogramLightResult = '否';
    histogramLightScore = 0;
    if lowLightPixelsRatio > histogramThreshold
        histogramLightResult = '是';
        histogramLightScore = (lowLightPixelsRatio - histogramThreshold) / histogramThreshold;
    end

    if size(image, 3) == 3
        img_gray = rgb2gray(image);
    else
        img_gray = image;
    end

    F = fft2(double(img_gray));
    F_shifted = fftshift(F);
    magnitude_spectrum = log(1 + abs(F_shifted));

    [rows, cols] = size(img_gray);
    center_row = floor(rows / 2);
    center_col = floor(cols / 2);
    radius = min(center_row, center_col) / 4;
    [X, Y] = meshgrid(1:cols, 1:rows);
    distance = sqrt((X - center_col).^2 + (Y - center_row).^2);

    high_freq = distance > radius;
    high_freq_energy = sum(sum(abs(F_shifted(high_freq)).^2));
    total_energy = sum(sum(abs(F_shifted).^2));
    high_freq_ratio = high_freq_energy / total_energy;
    blurResult = '否';
    blurScore = 0;
    if high_freq_ratio < blurThreshold
        blurResult = '是';
        blurScore = abs(high_freq_ratio - blurThreshold) / blurThreshold;
    end

    % if max_d < comprehensiveScore
    %     max_d = comprehensiveScore;
    % end
    comprehensiveScore = comprehensiveScore / 27.1241;
    % comprehensiveScore=comprehensiveScore/max(comprehensiveScore);
    LightScore=max(meanLightScore, histogramLightScore); 

  if comprehensiveScore > LightScore && comprehensiveScore > blurScore
    Degraded_Image_Classification = 'color_cast';
elseif LightScore > comprehensiveScore && LightScore > blurScore
    Degraded_Image_Classification = 'low_light';
elseif blurScore > comprehensiveScore && blurScore > LightScore
    Degraded_Image_Classification = 'blur';
end

    fprintf(fileID, '%s,%s,%s,%s,%s,%s,%.4f,%.4f,%.4f,%.4f,%.4f,%.4f,%s,%.4f,%s\n', file_name,meanDiffResult, varRatioResult, comprehensiveResult, meanLightResult, histogramLightResult, meanDiffScore, varRatioScore, comprehensiveScore, meanLightScore, histogramLightScore, LightScore,blurResult, blurScore,Degraded_Image_Classification);
end  
fclose(fileID);


