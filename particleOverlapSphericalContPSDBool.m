% Verify if there is particle overlap between the particle "part" and all
% the candidates in "collisionCandidates" inside the cell with positions x,y,z
function collision = particleOverlapSphericalContPSDBool(L, partDiameter, xPart, yPart, zPart,...
                                                         candidateDiameter, xCandidate, yCandidate, zCandidate)
    collision = false;

    dx = abs(xCandidate - xPart);
    dy = abs(yCandidate - yPart);
    dz = abs(zCandidate - zPart);

    % If particles are collision candidates but too distant, they may be overlapping through a boundary.
    dx = sign(L-2*dx) .* (dx - 0.5*L) + 0.5*L;
    dy = sign(L-2*dy) .* (dy - 0.5*L) + 0.5*L;
    dz = sign(L-2*dz) .* (dz - 0.5*L) + 0.5*L;

    dist = sqrt(power(dx,2) + power(dy,2) + power(dz,2));
        
    % If distance between centers < diameter, then there is a collision
    if any(dist <= (partDiameter + candidateDiameter)/2 & dist > 0, "all")
        collision = true;
    end
end