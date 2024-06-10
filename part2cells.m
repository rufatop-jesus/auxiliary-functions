function f = part2cells(nCells, cellSize,partDiameter, x, y, z)
    
    xLowerCell = ceil((x-partDiameter) / cellSize);
    yLowerCell = ceil((y-partDiameter) / cellSize);
    zLowerCell = ceil((z-partDiameter) / cellSize);

    xUpperCell = ceil((x+partDiameter) / cellSize);
    yUpperCell = ceil((y+partDiameter) / cellSize);
    zUpperCell = ceil((z+partDiameter) / cellSize);

    f = zeros((xUpperCell-xLowerCell) * (yUpperCell-yLowerCell) * (zUpperCell-zLowerCell), 3);  
    
    n = 1;
    for i = xLowerCell : xUpperCell
        for j = yLowerCell : yUpperCell
            for k = zLowerCell : zUpperCell

                if i > nCells
                    xCell = i - nCells;
                elseif i < 1
                    xCell = i + nCells;
                else
                    xCell = i;
                end

                if j > nCells
                    yCell = j - nCells;

                elseif j < 1
                    yCell = j + nCells;

                else
                    yCell = j;
                end

                if k > nCells
                    zCell = k - nCells;

                elseif k < 1
                    zCell = k + nCells;

                else
                    zCell = k;
                end

                f(n,1) = xCell;
                f(n,2) = yCell;
                f(n,3) = zCell;
                n = n+1;
            end
        end
    end
end