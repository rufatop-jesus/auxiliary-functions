% Calculate the frenquency of collisions between particle sizes

function collisionFreq = initCollisionFreq(diameter, stokesDist, X, Y, Z, L, PBIndex)

    sampleSize = 100;
    collisionFreq = zeros(size(diameter,1)); % array that registers the number of collisions between the different particle sizes

    % For each size index s
    for s = 1 : size(collisionFreq,1)
        s
        tic
        % Select sampleSize particles with size index s
        particlesS = find(PBIndex == s);
        firstPart = particlesS(1);
        lastPart = particlesS(end);
        particlesS = randi([firstPart, lastPart],1,sampleSize);
        
        % For each other size index n greater than or qual s
        for n = 1 : size(collisionFreq,2)       
            
            minStokesDist = min([stokesDist(s,1) stokesDist(n,1)]);
            diameterAvg = (diameter(s,1) + diameter(n,1))/2;
            particlesN = find(PBIndex == n); % select all particles with size index n
    
            for pS = particlesS
            
                % Calculate the x, y and z distances with neighbors with indexes
                % smaller than p, hence with diameters higher than or equal p has
                dx = abs(X(particlesN) - X(pS));
                dy = abs(Y(particlesN) - Y(pS));
                dz = abs(Z(particlesN) - Z(pS));
            
                % Correct diatnces considering periodic boundary condition
                dx = sign(L-2*dx) .* (dx - 0.5*L) + 0.5*L;
                dy = sign(L-2*dy) .* (dy - 0.5*L) + 0.5*L;
                dz = sign(L-2*dz) .* (dz - 0.5*L) + 0.5*L;
                
                % Calculate the distance
                distances = sqrt(power(dx,2) + power(dy,2) + power(dz,2));
        
                % Accumulate the possible collisions between particle sizes
                if minStokesDist > L/2
                    collisionFreq(s,n) = collisionFreq(s,n) + power(2*minStokesDist / L, 3) * nnz(distances < minStokesDist + diameterAvg & distances > 0);
                else
                    collisionFreq(s,n) = collisionFreq(s,n) + nnz(distances < minStokesDist + diameterAvg & distances > 0);
                end
            end        
        end
        toc
    end
    collisionFreq = collisionFreq/ sampleSize;
end