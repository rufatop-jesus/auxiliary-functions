% From the center coordinates (x y z), calculate which cells make the cube
% circunscribed about the particle with diamter = partDiameter.

function cells = initPart2Cells(nCells, cellSize, partDiameter, x, y, z)
    
    xLowerCell = ceil((x-partDiameter) / cellSize);
    xUpperCell = ceil((x+partDiameter) / cellSize);
    xCells = xLowerCell : xUpperCell; % cells that the particle occupies in the x direction
    nXCells = numel(xCells); % number of cells that the particle occupies in the x direction

    yLowerCell = ceil((y-partDiameter) / cellSize);
    yUpperCell = ceil((y+partDiameter) / cellSize);
    yCells = yLowerCell : yUpperCell; % cells that the particle occupies in the y direction
    nYCells = numel(yCells); % number of cells that the particle occupies in the y direction

    zLowerCell = ceil((z-partDiameter) / cellSize);
    zUpperCell = ceil((z+partDiameter) / cellSize);
    zCells = zLowerCell : zUpperCell; % cells that the particle occupies in the z direction
    nZCells = numel(zCells); % number of cells that the particle occupies in the z direction   

    xCells = repmat(xCells, [1, nYCells*nZCells]);
    yCells = repmat(repelem(yCells,nXCells), [1, nZCells]);
    zCells = repelem(zCells, nXCells*nYCells);

    cells = [xCells; yCells; zCells];
    
    % Apply periodic boundary condition for collision detection
    cells(cells > nCells) = cells(cells > nCells) - nCells;
    cells(cells < 1) = cells(cells < 1) + nCells;

    % cells = cells - nCells * ( sign(cells - nCells) + abs(sign(cells - nCells)) ) / 2; 
    % cells = cells - nCells * ( sign(cells - 1) - abs(sign(cells - 1)) ) / 2; 
    
    % Transform subscript to linear inices
    cells = (cells(3,:)-1) * nCells*nCells + (cells(2,:)-1) * nCells + cells(1,:);
    % cells = sub2ind([nCells,nCells,nCells], cells(1,:), cells(2,:), cells(3,:));  
end