% From the center coordinates (x y z), calculate which cells make the cube
% circunscribed about the particle with diamter = partDiameter.

function cells = initPart2Cells(nCells, cellSize, partDiameter, x, y, z)
    
    xLowerCell = ceil((x-partDiameter) / cellSize);
    xUpperCell = ceil((x+partDiameter) / cellSize);
    xCells = xLowerCell : xUpperCell;

    yLowerCell = ceil((y-partDiameter) / cellSize);
    yUpperCell = ceil((y+partDiameter) / cellSize);
    yCells = yLowerCell : yUpperCell;

    zLowerCell = ceil((z-partDiameter) / cellSize);
    zUpperCell = ceil((z+partDiameter) / cellSize);
    zCells = zLowerCell : zUpperCell;

    cells = combvec(xCells, yCells, zCells);
    
    % Apply periodic boundary condition for collision detection
    cells(cells > nCells) = cells(cells > nCells) - nCells;
    cells(cells < 1) = cells(cells < 1) + nCells;
    
    % Build meshgrid with the cells
    cells = sub2ind([nCells,nCells,nCells], cells(1,:), cells(2,:), cells(3,:));  
end