% From the current cells occupied by the particle and the test movement,
% calculate how much cells will the particle move on each cartesian
% direction and returns the cells it will occupy.

function f = testPart2cells(cells, nCells, cellSize, Xtest, Ytest, Ztest, X, Y, Z)
    f = cells;
    
    xCellMove = ceil(Xtest / cellSize) - ceil(X / cellSize);
    yCellMove = ceil(Ytest / cellSize) - ceil(Y / cellSize);
    zCellMove = ceil(Ztest / cellSize) - ceil(Z / cellSize);
    
    f(:,1) = f(:,1) + xCellMove;
    f(:,2) = f(:,2) + yCellMove;
    f(:,3) = f(:,3) + zCellMove;

    if any(f>nCells,'all')
        f(f>nCells) = f(f>nCells) - nCells;
    end

    if any(f<1,'all')
        f(f<1) = f(f<1) + nCells;
    end
end
