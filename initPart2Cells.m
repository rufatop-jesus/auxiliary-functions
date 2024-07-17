% From the center coordinates (x y z), calculate which cells make the cube
% circunscribed about the particle with diamter = partDiameter.

function cells = initPart2Cells(nCells, cellSize, partDiameter, x, y, z)
    
    xLowerCell = ceil((x-partDiameter) / cellSize);
    yLowerCell = ceil((y-partDiameter) / cellSize);
    zLowerCell = ceil((z-partDiameter) / cellSize);

    xUpperCell = ceil((x+partDiameter) / cellSize);
    yUpperCell = ceil((y+partDiameter) / cellSize);
    zUpperCell = ceil((z+partDiameter) / cellSize);

    cells = zeros((xUpperCell-xLowerCell+1) * (yUpperCell-yLowerCell+1) * (zUpperCell-zLowerCell+1), 3); 

    cellIndex = 1;
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

                cells(cellIndex,1) = xCell;
                cells(cellIndex,2) = yCell;
                cells(cellIndex,3) = zCell;
                cellIndex = cellIndex+1;
            end
        end
    end
end