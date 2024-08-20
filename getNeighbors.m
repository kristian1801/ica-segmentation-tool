function neighbors = getNeighbors(point, imgSize)
    % Define relative positions of 8-connected neighbors
    neighborOffsets = [-1, -1; -1, 0; -1, 1; 0, -1; 0, 1; 1, -1; 1, 0; 1, 1];
    
    % Calculate absolute positions of neighbors
    neighbors = bsxfun(@plus, neighborOffsets, point);
    
    % Remove neighbors that are outside the image boundaries
    neighbors(any(neighbors < 1, 2) | neighbors(:,1) > imgSize(1) | neighbors(:,2) > imgSize(2), :) = [];
end
