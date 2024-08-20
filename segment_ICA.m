function segment_ICA(filePath,pathAndNameToSave)
%%
% SEGMENT_ICA Segment vessels in TOF MR images using a 2D region growing algorithm.
%
% This function provides an interface to segment vessels in Time-of-Flight (TOF)
% Magnetic Resonance (MR) imaging. It utilizes a 2D region growth algorithm with
% dynamic intensity thresholding to accurately determine the vessel area. The
% user can interact via a GUI to save the segmentation results, which includes
% diameter estimates of the vessels, specifically the Internal Carotid Artery (ICA).
%
% Inputs:
%   filePath            - String specifying the full path to the .nii file containing
%                         the TOF MR image.
%   pathAndNameToSave   - String specifying the path and the subject-specific name
%                         under which the segmentation results (.csv file) will be saved.
%
% Outputs:
%   A .csv file is created at the specified location when the user presses the 'save'
%   button in the GUI. This file contains summary statistics of the vessel segmentation,
%   specifically the diameter estimates.
%
% Example:
%   segment_ICA('path/to/image.nii', 'path/to/save/results/mySegmentation.csv')
%
% Notes:
%   The raw TOF image is preprocessed before segmentation. For more details on the
%   preprocessing steps, refer to the documentation of preprocessImageData(imageData).
%
% See also:
%   preprocessImageData
%   getNeighbors
%   illustrateICA
%   regionGrowing2D
%
% Author:
%   Kristian Larsen, 15/04/2024

%% Preprocessing steps for the ToF image and data storage in struct
% Start by loading the niftiimage provided by the user
try
    % Attempt to load the NIfTI file using the normal method
    nii = load_nii(filePath);
    imageData = nii.img;
catch
    % If the normal loading method fails, try the alternative untouch loading method
    try
        fprintf('Normal loading failed, attempting untouch load...\n');
        nii = load_untouch_nii(filePath);
        imageData = nii.img;
        fprintf('Untouch loading succeeded.\n');
    catch ME
        % If both loading methods fail, handle the error
        fprintf('Error loading image data: %s\n', ME.message);
        return;
    end
end

% Extract pixel dimensions from the NIfTI header to convert pixels to millimeters
pixel_dims = nii.hdr.dime.pixdim(2:4);

% Preprocess the TOF (Time-of-Flight) image data to enhance the visibility of vessels
imageData = preprocessImageData(imageData);

% Determine a threshold value that is assumed to capture the highest blood signals
% This threshold is set to the 98.5th percentile of image data values, targeting the
% most intense (likely vessel) regions in the image
threshold = prctile(imageData(:), 98.5);

% Suppress all pixel intensities below the determined threshold, setting them to zero
% to focus on high-intensity vessel-related regions
imageData(imageData < threshold) = 0;

% Begin setting up the GUI for image exploration

% Retrieve the dimensions of the preprocessed image data
[rows, cols, slices] = size(imageData);

% Initialize a structure to store various image-related data
data = struct();
data.file = filePath; % Store the file path of the image

% Create a logical array to hold the binary masks for vessel segmentation per slice
data.masks = false(rows, cols, slices); % Initialize 3D mask array with false

% Initialize an array to store diameter estimates for each slice
data.diameters = zeros(1, slices);

% Track which slices are included in the analysis with a logical array
includedSlices = true(1, slices);

%% Set up GUI element

% Initialize the main figure for the GUI with specified properties
fig = figure('Name', '3D Image Explorer', 'NumberTitle', 'off', ...
    'Toolbar', 'none', 'MenuBar', 'none', ...
    'Position', [100, 100, 800, 600]); % Adjust [left, bottom, width, height]


% Define a custom function to execute when the GUI window is requested to close
set(fig, 'CloseRequestFcn', @closeFig);

% Create axes for image display
ax = axes('Parent', fig, 'Units', 'normalized', 'Position', [0.1, 0.2, 0.8, 0.7]);

% Display the first slice initially
imgHandle = imshow(imageData(:,:,1), []);

% Adjust the axes position slightly to leave space for the slider
ax.Position = [0.1, 0.2, 0.75, 0.7];

% Create a vertical slider for slice navigation on the right side of the GUI
slider = uicontrol('Parent', fig, 'Style', 'slider', 'Units', 'normalized', ...
    'Position', [0.9, 0.2, 0.05, 0.7], ...
    'Value', 1, 'Min', 1, 'Max', slices, ...
    'SliderStep', [1/(slices-1) , 10/(slices-1)], ...
    'Callback', @sliderCallback);

% Create a vertical button panel on the left side of the GUI
% Button for adding ROI manually
manualAddRoiButton = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Manual Add ROI', ...
    'Units', 'normalized', 'Position', [0.01, 0.85, 0.15, 0.05], ...
    'Callback', @manualAddRoiCallback);

% Button for adding ROI automatically
autoRoiButton = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Auto Delineate ROI', ...
    'Units', 'normalized', 'Position', [0.01, 0.78, 0.15, 0.05], ...
    'Callback', @autoRoiCallback);

% Button for adding ROI
addRoiButton = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Add ROI', ...
    'Units', 'normalized', 'Position', [0.01, 0.71, 0.15, 0.05], ...
    'Callback', @addRoiCallback);

% Button for deleting the current ROI
deleteRoiButton = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Delete ROI', ...
    'Units', 'normalized', 'Position', [0.01, 0.64, 0.15, 0.05], ...
    'Callback', @deleteRoiCallback);

% Button for 3D illustration
illustrate3DButton = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Illustrate 3D', ...
    'Units', 'normalized', 'Position', [0.01, 0.57, 0.15, 0.05], ...
    'Callback', @illustrate3DCallback);

% Button for saving the segmentation results to a CSV file
saveCsvButton = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Save CSV', ...
    'Units', 'normalized', 'Position', [0.01, 0.50, 0.15, 0.05], ...
    'Callback', @saveCsvCallback);

% Button for calculating diameters
calcDiameterButton = uicontrol('Parent', fig, 'Style', 'pushbutton', 'String', 'Calc Diameter', ...
    'Units', 'normalized', 'Position', [0.01, 0.43, 0.15, 0.05], ...
    'Callback', @calcDiameterCallback);

% Display for diameter statistics
diameterDisplay = uicontrol('Parent', fig, 'Style', 'text', ...
    'Units', 'normalized', 'Position', [0.01, 0.36, 0.15, 0.05], ...
    'String', 'Diameter Stats will be shown here');

% Toggle button for enabling and disabling zoom
zoomToggleButton = uicontrol('Parent', fig, 'Style', 'togglebutton', 'String', 'Zoom', ...
    'Units', 'normalized', 'Position', [0.01, 0.29, 0.15, 0.05], ...
    'Callback', @zoomToggleCallback);

%% Define callback functions for all buttons

% Define the close request function
    function closeFig(src, ~)
        % Cleanup or final processing before closing the GUI
        data.masks = data.masks(:,:,includedSlices);
        data.diameters = data.diameters(data.diameters > 1);
        % Close the figure
        delete(src);
    end

% Function to toggle zoom on and off
    function zoomToggleCallback(src, ~)
        if get(src, 'Value') % If the toggle button is pressed
            zoom on; % Enable zooming
        else
            zoom off; % Disable zooming
        end
    end

    function sliderCallback(src, ~)
        sliceNum = round(src.Value);
        if includedSlices(sliceNum)
            % Update the displayed image without clearing the axes
            set(imgHandle, 'CData', imageData(:,:,sliceNum));

            % Remove any existing ROI plots before overlaying new ones
            % This can be done by deleting the children of the axes that are line objects
            delete(findobj(ax, 'Type', 'line'));

            % Check and display ROI mask for the current slice
            if any(data.masks(:,:,sliceNum), 'all') % Check if there is any ROI for the current slice
                hold on; % Hold on to overlay the mask
                mask = data.masks(:,:,sliceNum);
                [B,~] = bwboundaries(mask, 'noholes');
                for k = 1:length(B)
                    boundary = B{k};
                    plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
                end
                hold off;
            end
            drawnow;
        else
            updateSliceDisplay(); % Skip discarded slices
        end
    end

    function addRoiCallback(~, ~)
        sliceNum = round(slider.Value);
        img = imageData(:,:,sliceNum);
        fprintf('%d is the current slice number in add roi\n', sliceNum);

        % Let user select a seed point
        [x, y] = ginput(1); % Gets seed point from user click
        seed = round([y, x]);

        % Perform region growing to get the ROI mask
        mask = regionGrowing2D(img, seed);

        % Calculate equivalent diameter
        diameter = calculateDiameterFromMask(mask, pixel_dims);

        % Update array with the estimated volume parameter
        data.diameters(sliceNum) = diameter;
        data.masks(:,:,sliceNum) = mask;

        % Display the ROI mask overlayed on the current slice
        hold on; % Hold on to overlay the mask
        [B,~] = bwboundaries(mask, 'noholes');
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
        end
        hold off;
    end

    function deleteRoiCallback(~, ~)
        sliceNum = round(slider.Value);
        data.masks(:,:,sliceNum) = false(rows,cols); % Clear the mask for the current slice
        data.diameters(sliceNum) = NaN;
        updateDisplayWithoutMask(sliceNum); % Update display to show image without the ROI
    end

% Function to update the image display based on included slices
    function updateDisplayWithoutMask(sliceNum)
        if isempty(sliceNum)
            % If no more slices are included, disable the slider and button
            set([slider, discardButton], 'Enable', 'off');
        else
            delete(findobj(ax, 'Type', 'line'));
            set(slider, 'Value', sliceNum);
            set(imgHandle, 'CData', imageData(:,:,sliceNum));
        end
        drawnow;
    end


    function manualAddRoiCallback(~, ~)
        sliceNum = round(slider.Value);

        % Activate the polygon drawing tool on the current axes
        h = drawpolygon('Color', 'r');

        % Wait for the user to finish drawing the polygon
        customWait(h);

        % Convert the polygon to a mask and update the data structure
        mask = createMask(h, imgHandle);
        data.masks(:,:,sliceNum) = mask;

        % Calculate and update the diameter based on the new ROI
        diameter = calculateDiameterFromMask(mask, pixel_dims);
        data.diameters(sliceNum) = diameter;

        % Update the display to show the new ROI
        updateDisplayWithMask(sliceNum);
    end

    function customWait(hROI)
        % Wait for the user to double-click the ROI to finalize it
        wait(hROI);
    end

    function diameter = calculateDiameterFromMask(mask, pixel_dims)
        areaPixels = sum(mask(:)); % Number of pixels in the ROI
        areaMM2 = areaPixels * pixel_dims(1) * pixel_dims(2); % Convert area to mm²
        diameter = 2 * sqrt(areaMM2 / pi); % Diameter of a circle with the same area
        return
    end

    function updateDisplayWithMask(sliceNum)
        set(imgHandle, 'CData', imageData(:,:,sliceNum));
        hold on;
        % Display the ROI mask for the current slice
        mask = data.masks(:,:,sliceNum);
        [B,~] = bwboundaries(mask, 'noholes');
        for k = 1:length(B)
            boundary = B{k};
            plot(boundary(:,2), boundary(:,1), 'r', 'LineWidth', 2);
        end
        hold off;
        drawnow;
    end

    function calcDiameterCallback(~, ~)
        % Assume includedSlices and data.diameters are well-aligned
        filteredDiameters = data.diameters(includedSlices);
        % Further filter diameters based on a specific condition
        validDiameters = filteredDiameters(filteredDiameters > 1);

        if isempty(validDiameters)
            statsStr = 'No valid diameters found.';
        else
            % Calculate summary statistics
            meanDiameter = mean(validDiameters);
            sdDiameter = std(validDiameters);
            medianDiameter = median(validDiameters);
            minDiameter = min(validDiameters);
            maxDiameter = max(validDiameters);

            % Format the statistics into a string for display
            statsStr = sprintf('Mean: %.2f, SD: %.2f, Median: %.2f, min: %.2f, max %.2f', ...
                meanDiameter, sdDiameter, medianDiameter, minDiameter, maxDiameter);
        end

        % Update the static text control to display the statistics
        set(diameterDisplay, 'String', statsStr);
    end


    function autoRoiCallback(~, ~)
        try
            % Let user select a seed point on the current slice.
            fprintf('Select a seed point on the current slice.\n');
            [x, y] = ginput(1); % Gets seed point from user click
            if isempty(x) || isempty(y)
                disp('ROI selection cancelled.');
                return;
            end
            seed = round([y, x]);
            previousSeed = seed; % Initialize previousSeed with the first seed.
            currentSlice = round(slider.Value);
            % Process each slice using the regionGrowing2D function.
            for sliceNum = currentSlice:slices
                img = imageData(:,:,sliceNum);
                mask = regionGrowing2D(img, seed);
                if ~any(mask(:))
                    disp(['Failed to find ROI in slice ', num2str(sliceNum), '. Using previous seed.']);
                    seed = previousSeed; % Fallback to the previous successful seed if current failed.
                else
                    % Update data structure if ROI is successfully found.
                    props = regionprops(mask, 'Centroid');
                    if ~isempty(props) % Ensure properties are calculated.
                        seed = round([props.Centroid(2), props.Centroid(1)]); % Update seed based on new ROI centroid.
                        previousSeed = seed; % Update previousSeed with the new successful seed.
                    end

                    % Calculate equivalent diameter.
                    diameter = calculateDiameterFromMask(mask, pixel_dims);

                    % Update array with the estimated diameter.
                    data.diameters(sliceNum) = diameter;
                    data.masks(:,:,sliceNum) = mask;
                end
            end

            % Update display for the current slice after processing all slices.
            updateDisplayWithMask(round(slider.Value));

            % Notify the user that automatic ROI delineation is complete.
            disp('Auto ROI delineation completed for all slices.');
        catch ME
            disp(['An error occurred during automatic ROI delineation: ', ME.message]);
        end
    end

    function illustrate3DCallback(~, ~)
        illustrateICA(data.masks);
    end

    function saveCsvCallback(~, ~)
        % Check if diameters array is not empty
        if isempty(data.diameters)
            disp('Diameter data is empty.');
            return;
        end

        % Assume includedSlices and data.diameters are well-aligned
        filteredDiameters = data.diameters(includedSlices);
        % Further filter diameters based on a specific condition
        validDiameters = filteredDiameters(filteredDiameters > 1);

        % Calculate summary statistics
        meanDiameter = mean(validDiameters);
        stdDiameter = std(validDiameters);
        medianDiameter = median(validDiameters);
        minDiameter = min(validDiameters);
        maxDiameter = max(validDiameters);
        NoOfSlices = length(validDiameters);

        % Display summary statistics
        disp(['Mean Diameter: ', num2str(meanDiameter)]);
        disp(['Standard Deviation of Diameter: ', num2str(stdDiameter)]);
        disp(['Median Diameter: ', num2str(medianDiameter)]);
        disp(['Minimum Diameter: ', num2str(minDiameter)]);
        disp(['Maximum Diameter: ', num2str(maxDiameter)]);

        % Prepare data for CSV in a wide format
        summaryStatsTable = table(meanDiameter, stdDiameter, medianDiameter, minDiameter, maxDiameter, NoOfSlices, ...
            'VariableNames', {'MeanDiameter', 'StdDiameter', 'MedianDiameter', 'MinDiameter', 'MaxDiameter', 'NoOfSlices'});

        % Define CSV file name
        csvFileName = [pathAndNameToSave '_diameter_statistics.csv'];

        % Write summary statistics to CSV file in wide format
        try
            writetable(summaryStatsTable, csvFileName);
            disp(['Summary statistics saved to ', csvFileName]);
        catch ME
            disp(['Failed to save summary statistics to CSV. Error: ', ME.message]);
        end
    end

end
