function [cell2parts, part2cells, PBIndex, X, Y, Z] = initPositionsRSA(nPart, partNumberDistribution, diameter, L)
    
    % Initialize arrays
    PBIndex = zeros(nPart, 1);
    X = zeros(nPart, 1);
    Y = zeros(nPart, 1);
    Z = zeros(nPart, 1);
    
    % Mesh of cubic cells to speed-up collision detection
    nCells = ceil(L * sqrt(3) / diameter(end,1));
    cellSize = L / nCells;
    cell2parts = cell(nCells,nCells,nCells); % list of particles within a cell   
    part2cells = cell(nPart,1); % list of cells a particles occupies

    totAdParticles = 0; % total of adsorbed particles
    isCandidate = zeros(1, nPart); % list with 0s and 1s. If the element (1, i) == 1, then the particle i is a candidate for collision

    % For every particle diameter do
    for i = size(partNumberDistribution, 1) : -1 : 1
        i
        collisions = 0;
        adParticles = 0; % total of adsorbed particles of the current size

        % While added particles < total of particles in the current size do
        while adParticles < partNumberDistribution(i,1)
            
            part = totAdParticles + adParticles + 1 % current particle under adsorption
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
            cells = initPart2cells(nCells, cellSize, diameter(i,1), xTest, yTest, zTest);
            
            % For every cell in cells find the particles that occupy it.
            for c = 1 : size(cells,1)
                isCandidate(cell2parts{cells(c,1), cells(c,2), cells(c,3)}) = 1;
            end
            collisionCandidates = find(isCandidate == 1); % list with the indexes of candidate particles for collision
            isCandidate(collisionCandidates) = 0; % clean the isCandidate array
            
            collisionList = particleOverlapSphericalContPSD(collisionCandidates, cellSize, part, PBIndex, diameter, L, X(1:part), Y(1:part), Z(1:part));

            % If there is no collision do
            if numel(collisionList) == 0
                % Register the cells the particle under adsorption occupy in the part2cells array
                part2cells{part,1} = cells; 
                
                % Register the particle under adsorption in the cell2parts array 
                for c = 1 : size(cells,1)
                    cell2parts{cells(c,1), cells(c,2), cells(c,3)} = [cell2parts{cells(c,1), cells(c,2), cells(c,3)} part];
                end

                adParticles = adParticles + 1; % sum one in the adsorbed particles counter

                
            else
                collisions = collisions + 1;
            end
            
        end
        % Volume percentage occupied by adsorbed particles
        dot(partNumberDistribution(end : -1 : i, 1) , (pi/6) * power(diameter(end : -1: i, 1),3)) / power(L,3) 
        % Collision rate
        adParticles / (collisions + adParticles);

        totAdParticles = totAdParticles + adParticles; % sum the counter of the size to the aggregate counter
    end
end