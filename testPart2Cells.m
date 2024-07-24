% From the current cells occupied by the particle and the test movement,
% calculate how much cells will the particle move on each cartesian
% direction and returns the cells it will occupy.

function newCells = testPart2Cells(cells, nCells, cellSize, XTest, YTest, ZTest, X, Y, Z)
    
    % Initialize the newCells arrays equal to the current cells
    newCells = cells;
    
    % Calculate how many cells were passed by the random movement in each
    % cartesian direction
    xCellMove = ceil(XTest / cellSize) - ceil(X / cellSize);
    yCellMove = ceil(YTest / cellSize) - ceil(Y / cellSize);
    zCellMove = ceil(ZTest / cellSize) - ceil(Z / cellSize);
    
    % Sum(subtract) the number of passed by cells to the current cells to reach the
    % new cells
    newCells(:,1) = newCells(:,1) + xCellMove;
    newCells(:,2) = newCells(:,2) + yCellMove;
    newCells(:,3) = newCells(:,3) + zCellMove;
    
    % Correct for cells outside the upper boundary with a periodic boundary
    % condition
    newCells(newCells>nCells) = newCells(newCells>nCells) - nCells;
    
    % Correct for cells outside the lower boundary with a periodic boundary
    % condition
    newCells(newCells<1) = newCells(newCells<1) + nCells;
end