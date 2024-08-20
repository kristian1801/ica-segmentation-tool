function mask = regionGrowing2D(img, seed, maxRegionSize)
% REGIONGROWING2D Performs region growing segmentation on a 2D image.
%
% This function applies a region growing algorithm to segment an area in a 2D image
% starting from a user-defined seed point. The segmentation grows until the region
% reaches a user-specified size or the intensity difference threshold is exceeded.
% The threshold is adaptive, depending on the local mean and standard deviation.
%
% Inputs:
%   img           - A 2D grayscale image matrix.
%   seed          - A 1x2 vector [row, col] indicating the seed point for the region
%                   growing algorithm.
%   maxRegionSize - An optional parameter specifying the maximum number of pixels
%                   in the region. 
%
% Outputs:
%   mask          - A logical mask of the same size as 'img', where pixels that are
%                   part of the region are 1 (true) and others are 0 (false). The mask
%                   is empty if the region grows beyond the maximum allowed size.
%
% Example:
%   img = imread('sample.jpg');
%   seed = [100, 150];
%   mask = regionGrowing2D(img, seed);
%
% Notes:
%   The region growing criteria is based on the intensity difference to the seed
%   value not exceeding the adaptive threshold, which is calculated as:
%   `localMean - 0.5 * localStd`.
%
%   The function stops and returns an empty mask if the region size exceeds
%   'maxRegionSize', indicating the region of interest might not be accurately
%   identified.
%
% Author:
%   Kristian Larsen, 15/04/2024

% Validate input arguments
if nargin < 3
    maxRegionSize = numel(img) / 90;
end

% Initialize the region mask and other variables
mask = false(size(img));
seedValue = img(seed(1), seed(2));

% List to manage the growth points
pointsToCheck = {seed};
regionSize = 0; % Initialize region size counter

% Region growing algorithm
while ~isempty(pointsToCheck) && regionSize <= maxRegionSize
    currentPoint = pointsToCheck{1};
    pointsToCheck(1) = [];
    
    % Retrieve neighbors of the current point
    neighbors = getNeighbors(currentPoint, size(img));
    
    % Process each neighbor
    for i = 1:size(neighbors, 1)
        neighborPoint = neighbors(i, :);
        x = neighborPoint(1);
        y = neighborPoint(2);
        
        % Check if the neighbor can be added to the region
        if ~mask(x, y)
            localMean = mean2(img(max(x-1,1):min(x+1,end), max(y-1,1):min(y+1,end)));
            localStd = std2(img(max(x-1,1):min(x+1,end), max(y-1,1):min(y+1,end)));
            adaptiveThreshold = localMean - localStd; % Adaptive formula
            
            if abs(img(x, y) - seedValue) <= adaptiveThreshold
                mask(x, y) = true;
                pointsToCheck{end+1} = neighborPoint;
                regionSize = regionSize + 1; % Increment the region size counter
            end
        end
    end
end

% Handle the case where the region size exceeds the limit
if regionSize > maxRegionSize
    mask(:) = false;
    disp('Maximum region size exceeded, invalidating ROI.');
end

end
