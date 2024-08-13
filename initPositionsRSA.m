function [cell2parts, nCells, cellSize, PBIndex, X, Y, Z] = initPositionsRSA(nPart, partNumberDistribution, diameter, avgDiameter, L)
    
    % Initialize arrays
    PBIndex = zeros(nPart, 1);
    X = zeros(nPart, 1);
    Y = zeros(nPart, 1);
    Z = zeros(nPart, 1);

    % Mesh of cubic cells to speed-up collision detection
    nCells = ceil(L / avgDiameter);
    cellSize = L / nCells;
    cell2parts = cell(nCells,nCells,nCells); % list of particles within a cell

    totAdParticles = 0; % total of adsorbed particles

    % For every particle diameter do
    for i = size(partNumberDistribution, 1) : -1 : 1
        tic
        i
        adParticles = 0; % total of adsorbed particles of the current size
        partDiameter = diameter(i,1);

        % Minimum distance to even be considered a collision candidate
        collisionDist = (partDiameter + diameter(:,1)) / 2;

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

            collisionCandidates = horzcat(cell2parts{cells(1,:)});

            dx = abs(X(collisionCandidates) - X(part));
            dy = abs(Y(collisionCandidates) - Y(part));
            dz = abs(Z(collisionCandidates) - Z(part));

            dx = sign(L-2*dx) .* (dx - 0.5*L) + 0.5*L;
            dy = sign(L-2*dy) .* (dy - 0.5*L) + 0.5*L;
            dz = sign(L-2*dz) .* (dz - 0.5*L) + 0.5*L;

            dist = sqrt(power(dx,2) + power(dy,2) + power(dz,2));

            if any(dist <= collisionDist(PBIndex(collisionCandidates)), "all")
                continue
            
            else
                % Register the particle under adsorption in the cell2parts array 
                for c = 1 : numel(cells)
                    cell2parts{cells(1,c)} = [cell2parts{cells(1,c)} part];
                end

                adParticles = adParticles + 1; % sum one in the adsorbed particles counter       
            end
         
        end
        totAdParticles = totAdParticles + adParticles; % sum the counter of the size to the aggregate counter
        toc
    end
    
end