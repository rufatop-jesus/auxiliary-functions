function [cell2parts, part2cells, nCells, cellSize, PBIndex, X, Y, Z] = initPositionsRSA(nPart, partNumberDistribution, diameter, L)
    
    % Initialize arrays
    PBIndex = zeros(nPart, 1);
    X = zeros(nPart, 1);
    Y = zeros(nPart, 1);
    Z = zeros(nPart, 1);
    
    % PND weighted average diameter
    % avgDiameter = dot(partNumberDistribution,diameter) / nPart;
    avgDiameter = max(diameter,[],"all");

    % Mesh of cubic cells to speed-up collision detection
    nCells = ceil(L * sqrt(3) / avgDiameter);
    cellSize = L / nCells;
    cell2parts = zeros(power(nCells,3), ceil(nPart/power(nCells,3))); % list of particles within a cell. nPart/power(nCells,3) is the estimated max number of particles within a cell
    nAddedParts = zeros(size(cell2parts,1),1); % array with the number of added particles on each cell of the cell2parts array.

    part2cells = zeros(nPart, ceil(power(max(diameter,[],"all"),3)/ power(cellSize,3))); % list of cells a particles occupies. power(max(diameter,[],"all"),3)/ power(cellSize,3) is the estmated max number of cells a particle can occupy

    totAdParticles = 0; % total of adsorbed particles

    % For every particle diameter do
    for i = size(partNumberDistribution, 1) : -1 : 1
        i
        adParticles = 0; % total of adsorbed particles of the current size

        % While added particles < total of particles in the current size do
        while adParticles < partNumberDistribution(i,1)
            
            part = totAdParticles + adParticles + 1; % current particle under adsorption
            PBIndex(part) = i; % index of this particle in the Particle Book

            % Select a random cell within the domain for the particle under adsorption
            xCellTest = randi(nCells);
            yCellTest = randi(nCells);
            zCellTest = randi(nCells);

            % Generate a random position within the domain for the particle under adsorption
            xTest = (xCellTest - rand) * cellSize;
            yTest = (yCellTest - rand) * cellSize;
            zTest = (zCellTest - rand) * cellSize;

            % Register the random position in the respective X,Y and Z arrays
            X(part) = xTest;
            Y(part) = yTest;
            Z(part) = zTest;
            
            % Identify the cells occupied by the particle under adsorption
            cells = initPart2Cells(nCells, cellSize, diameter(i,1), xTest, yTest, zTest);
            cells = sub2ind([nCells,nCells,nCells], cells(:,1), cells(:,2), cells(:,3));
            
            % collisionCandidates = [];
            % % For every cell in cells find the particles that occupy it.
            % for testingCell = 1 : size(cellsInd,1)
            %     if nAddedParts(cellsInd(testingCell)) > 0
            %         collisionCandidates = [collisionCandidates cell2parts(cellsInd(testingCell), 1:nAddedParts(cellsInd(testingCell)))];
            %     end
            % end
            collisionCandidates = unique(nonzeros(cell2parts(cells,:)));

            % Test for collision
            collision = particleOverlapSphericalContPSDBool(L, diameter(PBIndex(part),1), X(part), Y(part), Z(part),...
                                                            diameter(PBIndex(collisionCandidates),1),...
                                                            X(collisionCandidates), Y(collisionCandidates), Z(collisionCandidates));
            
            % If there is no collision do
            if not(collision)

                % Register the cells the particle under adsorption occupy in the part2cells array
                part2cells(part, 1 : size(cells,1)) = cells'; 

                % Register the particle under adsorption in the cell2parts array 
                for adsorbingCell = 1 : size(cells,1)
                    cell2parts(cells(adsorbingCell,1), nAddedParts(adsorbingCell)+1) = part;
                end

                adParticles = adParticles + 1; % sum one in the adsorbed particles counter       
                nAddedParts(cells,1) = nAddedParts(cells,1) + 1; % sum one on the cells the particle 'part' occupies
            end
         
        end
        totAdParticles = totAdParticles + adParticles; % sum the counter of the size to the aggregate counter
    end
end