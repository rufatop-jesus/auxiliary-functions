% Verify if there is particle overlap between the particle "part" and all
% the candidates in "collisionCandidates" inside the cell with positions x,y,z
function collisionList = particleOverlapSphericalContPSDList(collisionCandidates, L, partDiameter, xPart, yPart, zPart,...
                                                            candidateDiameter, xCandidate, yCandidate, zCandidate)
        
    dx = abs(xCandidate - xPart);
    dy = abs(yCandidate - yPart);
    dz = abs(zCandidate - zPart);

    % If particles are collision candidates but too distant, they may be overlapping through a boundary.
    dx = sign(L-2*dx) .* (dx - 0.5*L) + 0.5*L;
    dy = sign(L-2*dy) .* (dy - 0.5*L) + 0.5*L;
    dz = sign(L-2*dz) .* (dz - 0.5*L) + 0.5*L;
    
    % Minimum distance to even be considered a collision candidate
    collisionDist = (partDiameter + candidateDiameter) / 2;
    
    dist = sqrt(power(dx,2) + power(dy,2) + power(dz,2));

    actualCollisions = (dist < collisionDist & dist ~= 0);
    
    collisionList = collisionCandidates(actualCollisions);
end