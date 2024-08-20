function imageData = preprocessImageData(imageData)
% PREPROCESSIMAGEDATA Enhance the visibility of vasculature in 3D medical images.
%
% This function preprocesses 3D medical image data to enhance the visibility of
% vasculature structures by applying several image processing techniques. These
% include background subtraction, anisotropic diffusion, and a top-hat transform,
% aimed at improving the contrast and clarity of vascular features in the image.
%
% Input:
%   imageData - A 3D matrix representing the grayscale intensity values of a
%               medical image (e.g., MRI, CT). It must be a 3D array.
%
% Output:
%   imageData - The processed 3D image matrix with enhanced vasculature visibility.
%
% Usage:
%   processedData = preprocessImageData(rawImageData);
%
% Notes:
%   The function throws an error if the input is not a 3D array. This is crucial to
%   ensure the data integrity for the specific processing techniques used, which
%   are tailored for 3D volumetric data.
%
% Example of preprocessing steps:
%   1. Background subtraction using morphological opening with a disk-shaped
%      structuring element of radius 15 to eliminate uneven background illumination.
%   2. Anisotropic diffusion filtering to reduce image noise while preserving edge
%      sharpness, enhancing the visibility of smaller vessels.
%   3. Top-hat filtering using the same structuring element to highlight small
%      structures darker than their surroundings.
%
% Author:
%   Kristian Larsen, 15/04/2024

if ndims(imageData) ~= 3
    error('Input must be a 3D image.');
end


% Background subtraction to remove low-frequency background variations
background = imopen(imageData, strel('disk', 15));
imageData = imsubtract(imageData, background);
% Apply anisotropic diffusion filter to reduce noise and preserve edges
imageData = imdiffusefilt(imageData);
% Use top-hat transform to enhance brightness of small elements
se = strel('disk', 15);
imageData = imtophat(imageData, se);
end
