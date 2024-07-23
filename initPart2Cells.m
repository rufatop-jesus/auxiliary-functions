% From the center coordinates (x y z), calculate which cells make the cube
% circunscribed about the particle with diamter = partDiameter.

function cells = initPart2Cells(nCells, cellSize, partDiameter, x, y, z)
    
    xLowerCell = ceil((x-partDiameter) / cellSize);
    yLowerCell = ceil((y-partDiameter) / cellSize);
    zLowerCell = ceil((z-partDiameter) / cellSize);

    xUpperCell = ceil((x+partDiameter) / cellSize);
    yUpperCell = ceil((y+partDiameter) / cellSize);
    zUpperCell = ceil((z+partDiameter) / cellSize);
    
    % Build meshgrid with the cells
    [xCells, yCells, zCells] = meshgrid(xLowerCell:xUpperCell, yLowerCell:yUpperCell, zLowerCell:zUpperCell);
    
    % Reshape the meshgrid to a n-by-3 matrix
    cells = [reshape(xCells,[],1), reshape(yCells,[],1), reshape(zCells,[],1)];
    
    % Applt periodic boundary condition for collision detection
    cells(cells > nCells) = cells(cells > nCells) - nCells;
    cells(cells < 1) = cells(cells < 1) + nCells;
end