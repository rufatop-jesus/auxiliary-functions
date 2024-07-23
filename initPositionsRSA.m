        function [cell2parts, part2cells, PBIndex, X, Y, Z] = initPositionsRSA(nPart, partNumberDistribution, diameter, L)
    
    % Initialize arrays
    PBIndex = zeros(nPart, 1);
    X = zeros(nPart, 1);
    Y = zeros(nPart, 1);
    Z = zeros(nPart, 1);
    
    % PND weighted average diameter
    avgDiameter = dot(partNumberDistribution,diameter) / nPart;
    % avgDiameter = max(diameter,[],"all");

    % Mesh of cubic cells to speed-up collision detection
    nCells = ceil(L * sqrt(3) / avgDiameter);
    cellSize = L / nCells;
    cell2parts = cell(nCells,nCells,nCells); % list of particles within a cell   
    part2cells = cell(nPart,1); % list of cells a particles occupies

    totAdParticles = 0; % total of adsorbed particles

    % For every particle diameter do
    for i = size(partNumberDistribution, 1) : -1 : 1
        i
        adParticles = 0; % total of adsorbed particles of the current size
        collision = 0;

        % While added particles < total of particles in the current size do
        while adParticles < partNumberDistribution(i,1)
            
            part = totAdParticles + adParticles + 1; % current particle under adsorption
            PBIndex(part) = i; % index of this particle in the Particle Book

            % Select a random cell within the domain for the particle under adsorption
            xCellTest = randi(nCells);
            yCellTest = randi(nCells);
            zCellTest = randi(nCells);

            % while numel(cell2parts{xCellTest, yCellTest, zCellTest}) > 0
            %     xCellTest = randi(nCells);
            %     yCellTest = randi(nCells);
            %     zCellTest = randi(nCells);
            % end

            % Generate a random position within the domain for the particle under adsorption
            xTest = (xCellTest - rand) * cellSize;
            yTest = (yCellTest - rand) * cellSize;
            zTest = (zCellTest - rand) * cellSize;

            % Register the random position in the respective X,Y and Z arrays
            X(part) = xTest;
            Y(part) = yTest;
            Z(part) = zTest;
            
            % Identify the cells occupied by the particle under adsorption
            cellsSub = initPart2Cells(nCells, cellSize, diameter(i,1), xTest, yTest, zTest);
            cellsInd = sub2ind(size(cell2parts), cellsSub(:,1), cellsSub(:,2), cellsSub(:,3));
            
            % For every cell in cells find the particles that occupy it.
            collisionCandidates = cell2parts{cellsInd};
            % collisionCandidates = unique(collisionCandidates);

            isCollision = particleOverlapSphericalContPSD(L, diameter(PBIndex(part),1), X(part), Y(part), Z(part),...
                                                          diameter(PBIndex(collisionCandidates),1), X(collisionCandidates), Y(collisionCandidates), Z(collisionCandidates));
            
            % If there is no collision do
            if not(isCollision)
                
                % Register the cells the particle under adsorption occupy in the part2cells array
                part2cells{part,1} = cellsSub; 

                % % Register the particle under adsorption in the cell2parts array 
                for adsorbingCell = 1 : size(cellsSub,1)
                    cell2parts{cellsSub(adsorbingCell,1), cellsSub(adsorbingCell,2), cellsSub(adsorbingCell,3)} = [cell2parts{cellsSub(adsorbingCell,1), cellsSub(adsorbingCell,2), cellsSub(adsorbingCell,3)} part];
                end

                adParticles = adParticles + 1; % sum one in the adsorbed particles counter                  
            else
                collision = collision + 1;
            end
         
        end
        collision / adParticles
        totAdParticles = totAdParticles + adParticles; % sum the counter of the size to the aggregate counter
    end
end