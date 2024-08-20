function illustrateICA(ICA)
% ILLUSTRATEICA Visualizes the 3D structure of the Internal Carotid Artery (ICA) using isosurface plotting.
%
% This function generates a 3D visualization of the Internal Carotid Artery from a 
% 3D MRI data set. It uses an isosurface plot to represent the structure of the ICA 
% based on a specified intensity threshold. The plot is rendered in red without edge 
% colors to emphasize the surface.
%
% Input:
%   ICA - A 3D matrix representing the grayscale intensity values of a medical
%         image, specifically targeting the region containing the ICA.
%
% Usage:
%   illustrateICA(ICAdata);
%
% Notes:
%   The input matrix must be a 3D array, otherwise, the function throws an error.
%   Adjust the view and lighting settings as necessary to best visualize the anatomy.
%
% Example:
%   % Assuming 'ICAdata' is a preloaded 3D matrix of an ICA scan
%   illustrateICA(ICAdata);
%
% This function sets up the figure, configures the viewing angle, extracts the
% isosurface for the specified intensity threshold, and applies lighting to enhance
% visibility. Labels are added to each axis for clarity.
%
% See also:
%   isosurface, patch, light, lighting

    if ndims(ICA) ~= 3
        error('Input must be a 3D image.');
    end

    % Create a new figure window with high resolution
    fig = figure('Color', 'k', 'Position', [100 100 1200 900]);

    % Extract the isosurface data using multiple threshold values
    thresholds = [0.1, 0.3, 0.5];
    colors = [1 0.2 0.2; 0.8 0 0; 0.6 0 0];  % Different shades of red

    hold on;
    for i = 1:length(thresholds)
        p = patch(isosurface(ICA, thresholds(i)));
        isonormals(ICA, p);
        p.FaceColor = colors(i,:);
        p.EdgeColor = 'none';
        p.FaceAlpha = 0.8 - (i-1)*0.2;  % Decreasing opacity for inner layers
        p.SpecularStrength = 0.7;
        p.SpecularExponent = 5;
    end

    % Enhance lighting
    camlight('headlight');
    camlight('right');
    lighting gouraud;

    % Add a subtle ambient light
    light('Style', 'infinite', 'Position', [-1 -1 -1], 'Color', [0.3 0.3 0.4]);

    % Improve the appearance
    axis tight;
    axis off;
    view(45, 30);  % Adjust view angle for best presentation

    % Set background to gradient
    set(gca, 'Color', 'none');
    colormap(flipud(gray));  % For background gradient
    
    % Create gradient background
    hold on;
    fill3([0 1 1 0], [0 0 1 1], [-1 -1 -1 -1], [0 0 1 1], 'EdgeColor', 'none', 'FaceColor', 'interp');
    
    % Adjust camera properties for depth
    camva(8);  % Adjust camera view angle
    material dull;

    % Enhance overall lighting
    h = light('Position', [-1 0 0], 'Style', 'infinite');
    lightangle(-45,30);
    set(h, 'Color', [.6 .2 .2]);

    % Save high-resolution image
    %exportgraphics(fig, outputFileName, 'Resolution', 300);
    
    %close(fig);  % Close the figure after saving

end
