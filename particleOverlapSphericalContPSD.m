% Verify if there is particle overlap between the particle "part" and all
% the candidates in "collisionCandidates" inside the cell with positions x,y,z
function collision = particleOverlapSphericalContPSD(L, partDiameter, xPart, yPart, zPart,...
                                                    candidateDiameter, xCandidate, yCandidate, zCandidate)
    collision = false;
    
    % For all the collision candidates inside the cell do
    for i = 1 : size(candidateDiameter,1)
        
        dx = abs(xCandidate(i) - xPart);
        dy = abs(yCandidate(i) - yPart);
        dz = abs(zCandidate(i) - zPart);
 
        % If particles are in the same cell but too distant, they may be overlapping through a boundary.
        if dx > L-dx
            dx = L-dx;
        end

        if dy > L-dy
            dy = L-dy;
        end

        if dz > L-dz
            dz = L-dz;
        end
        
        dist = sqrt(power(dx,2) + power(dy,2) + power(dz,2)); % distance between centers
        
        % If distance between centers < diameter, then there is a collision
        if and(dist < (partDiameter + candidateDiameter(i))/2, dist > 0)
            collision = true;
            return
        end
    end
end