clc;    % Clear the command window.
close all;  % Close all figures (except those of imtool.)
imtool close all;  % Close all imtool figures if you have the Image Processing Toolbox.
clear;  % Erase all existing variables. Or clearvars if you want.
workspace;  % Make sure the workspace panel is showing.
format short g;
format compact;
folder=pwd;
fontSize = 20;
baseFileName = 'ag-no7_m008.tif';
% Get the full filename, with path prepended.
fullFileName = fullfile(folder, baseFileName);
% Check if file exists.
if ~exist(fullFileName, 'file')
  % File doesn't exist -- didn't find it there.  Check the search path for it.
  fullFileNameOnSearchPath = baseFileName; % No path this time.
  if ~exist(fullFileNameOnSearchPath, 'file')
    % Still didn't find it.  Alert user.
    errorMessage = sprintf('Error: %s does not exist in the search path folders.', fullFileName);
    uiwait(warndlg(errorMessage));
    return;
  end
end
grayImage = imread(fullFileName);
% Save this figure handle.
hFig1 = gcf;
% Get the dimensions of the image.  
% numberOfColorBands should be = 1.
[rows, columns, numberOfColorBands] = size(grayImage);
if numberOfColorBands > 1
  % It's not really gray scale like we expected - it's color.
  % Convert it to gray scale by taking only the green channel.
  grayImage = grayImage(:, :, 2); % Take green channel.
end
% Display the original gray scale image.
subplot(1, 3, 1);
imshow(grayImage, []);
title('Original Grayscale Image', 'FontSize', fontSize);
% Enlarge figure to full screen.
set(gcf, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Give a name to the title bar.
set(gcf, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off') 
binaryImage = grayImage < 128;
% Display the binary image.
subplot(1, 3, 2);
imshow(binaryImage, []);
title('Binary Image', 'FontSize', fontSize);
% Display the binary image again where we can put boxes over it.
subplot(1, 3, 3);
imshow(binaryImage, []);
hold on;
title('Binary Image', 'FontSize', fontSize);
% Create a figure for the cropped images.
hCropped = figure;
% Do connected components analysis on it.
cc = bwconncomp(binaryImage);
% Measure the bounding box of all blobs.
measurements = regionprops(cc, 'BoundingBox');
fprintf('Found %d regions\n', cc.NumObjects);
numSkinnyRegions = 0;
for k = 1 : cc.NumObjects
  figure(hFig1); % Switch to figure 1.
  thisBB = measurements(k).BoundingBox
  % Draw a box around the region in cyan.
  hRect = rectangle('Position', thisBB, 'EdgeColor', 'c', 'LineWidth', 3);
  aspectRatio(k) = thisBB(4)/thisBB(3);
  if (thisBB(4) <= 3 || thisBB(3) <= 3) && (aspectRatio(k) > 4 || aspectRatio(k) < 1/4)
    numSkinnyRegions = numSkinnyRegions + 1;
    % Save it to a cell array, just in case we want to use it after the loop is done.
    croppedImages{numSkinnyRegions} = imcrop(binaryImage, thisBB);    
    % Draw skinny regions in a different color
    delete(hRect); % Get rid of old one.
    hRect = rectangle('Position', thisBB, 'EdgeColor', 'r', 'LineWidth', 3);  
    % Switch to figure 2
    figure(hCropped);
    subplot(2, 3, numSkinnyRegions);
    imshow(croppedImages{numSkinnyRegions}, []);
    caption = sprintf('Blob #%d', numSkinnyRegions);
    title(caption, 'FontSize', fontSize);    
  end
end
% Enlarge figure to full screen.
set(hCropped, 'Units', 'Normalized', 'OuterPosition', [0 0 1 1]);
% Give a name to the title bar.
set(hCropped, 'Name', 'Demo by ImageAnalyst', 'NumberTitle', 'Off')