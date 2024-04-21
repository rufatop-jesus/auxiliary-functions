% Identify candidates to apply the particle overlap detection.
function [cd, xCd, yCd, zCd, alphaCd, betaCd, gammaCd] = identifyCandidates(xPOverlap,xNOverlap,yPOverlap,yNOverlap,zPOverlap,zNOverlap,...
                                                                            part,maxR,Lx,Ly,Lz,X,Y,Z,alpha,beta,gamma)
    % Identifying candidates within the domain
    cd = find(sqrt(power(X-X(part), 2) + power(Y-Y(part), 2) + power(Z-Z(part), 2)) < 2*maxR);
    xCd = X(cd);
    yCd = Y(cd);
    zCd = Z(cd);
    
    % Identifying candidates at the periodic boundaries on the faces, edges and vertices. XP = positive x,
    % XN = negative x and the same applies to y and z.
    if xPOverlap
        newCandidates = find(sqrt(power(X+Lx-X(part), 2) + power(Y-Y(part), 2) + power(Z-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) + Lx];
        yCd = [yCd ; Y(newCandidates)];
        zCd = [zCd ; Z(newCandidates)];
    end
    
    if xNOverlap
        newCandidates = find(sqrt(power(X-Lx-X(part), 2) + power(Y-Y(part), 2) + power(Z-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) - Lx];
        yCd = [yCd ; Y(newCandidates)];
        zCd = [zCd ; Z(newCandidates)];
    end

    if yPOverlap
        newCandidates = find(sqrt(power(X-X(part), 2) + power(Y+Ly-Y(part), 2) + power(Z-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates)];
        yCd = [yCd ; Y(newCandidates) + Ly];
        zCd = [zCd ; Z(newCandidates)];
    end
    
    if yNOverlap
        newCandidates = find(sqrt(power(X-X(part), 2) + power(Y-Ly-Y(part), 2) + power(Z-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates)];
        yCd = [yCd ; Y(newCandidates) - Ly];
        zCd = [zCd ; Z(newCandidates)];
    end

    if zPOverlap
        newCandidates = find(sqrt(power(X-X(part), 2) + power(Y-Y(part), 2) + power(Z+Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates)];
        yCd = [yCd ; Y(newCandidates)];
        zCd = [zCd ; Z(newCandidates) + Lz];
    end
    
    if zNOverlap
        newCandidates = find(sqrt(power(X-X(part), 2) + power(Y-Y(part), 2) + power(Z-Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates)];
        yCd = [yCd ; Y(newCandidates)];
        zCd = [zCd ; Z(newCandidates) - Lz];
    end

    if xPOverlap && yPOverlap
        newCandidates = find(sqrt(power(X+Lx-X(part), 2) + power(Y+Ly-Y(part), 2) + power(Z-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) + Lx];
        yCd = [yCd ; Y(newCandidates) + Ly];
        zCd = [zCd ; Z(newCandidates)];
    end

    if xNOverlap && yPOverlap
        newCandidates = find(sqrt(power(X-Lx-X(part), 2) + power(Y+Ly-Y(part), 2) + power(Z-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) - Lx];
        yCd = [yCd ; Y(newCandidates) + Ly];
        zCd = [zCd ; Z(newCandidates)];
    end

    if xPOverlap && yNOverlap
        newCandidates = find(sqrt(power(X+Lx-X(part), 2) + power(Y-Ly-Y(part), 2) + power(Z-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) + Lx];
        yCd = [yCd ; Y(newCandidates) - Ly];
        zCd = [zCd ; Z(newCandidates)];
    end
    
    if xNOverlap && yNOverlap
        newCandidates = find(sqrt(power(X-Lx-X(part), 2) + power(Y-Ly-Y(part), 2) + power(Z-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) - Lx];
        yCd = [yCd ; Y(newCandidates) - Ly];
        zCd = [zCd ; Z(newCandidates)];
    end

    if xPOverlap && zPOverlap
        newCandidates = find(sqrt(power(X+Lx-X(part), 2) + power(Y-Y(part), 2) + power(Z+Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) + Lx];
        yCd = [yCd ; Y(newCandidates)];
        zCd = [zCd ; Z(newCandidates) + Lz];
    end

    if xNOverlap && zPOverlap
        newCandidates = find(sqrt(power(X-Lx-X(part), 2) + power(Y-Y(part), 2) + power(Z+Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) - Lx];
        yCd = [yCd ; Y(newCandidates)];
        zCd = [zCd ; Z(newCandidates) + Lz];
    end

    if xPOverlap && zNOverlap
        newCandidates = find(sqrt(power(X+Lx-X(part), 2) + power(Y-Y(part), 2) + power(Z-Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) + Lx];
        yCd = [yCd ; Y(newCandidates)];
        zCd = [zCd ; Z(newCandidates) - Lz];
    end
    
    if xNOverlap && zNOverlap
        newCandidates = find(sqrt(power(X-Lx-X(part), 2) + power(Y-Y(part), 2) + power(Z-Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) - Lx];
        yCd = [yCd ; Y(newCandidates)];
        zCd = [zCd ; Z(newCandidates) - Lz];
    end
    
    if yPOverlap && zPOverlap
        newCandidates = find(sqrt(power(X-X(part), 2) + power(Y+Ly-Y(part), 2) + power(Z+Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates)];
        yCd = [yCd ; Y(newCandidates) + Ly];
        zCd = [zCd ; Z(newCandidates) + Lz];
    end

    if yNOverlap && zPOverlap
        newCandidates = find(sqrt(power(X-X(part), 2) + power(Y-Ly-Y(part), 2) + power(Z+Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates)];
        yCd = [yCd ; Y(newCandidates) - Ly];
        zCd = [zCd ; Z(newCandidates) + Lz];
    end

    if yPOverlap && zNOverlap
        newCandidates = find(sqrt(power(X-X(part), 2) + power(Y+Ly-Y(part), 2) + power(Z-Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates)];
        yCd = [yCd ; Y(newCandidates) + Ly];
        zCd = [zCd ; Z(newCandidates) - Lz];
    end
    
    if yNOverlap && zNOverlap
        newCandidates = find(sqrt(power(X-X(part), 2) + power(Y-Ly-Y(part), 2) + power(Z-Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates)];
        yCd = [yCd ; Y(newCandidates) - Ly];
        zCd = [zCd ; Z(newCandidates) - Lz];
    end
    
    if xPOverlap && yPOverlap && zPOverlap
        newCandidates = find(sqrt(power(X+Lx-X(part), 2) + power(Y+Ly-Y(part), 2) + power(Z+Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) + Lx];
        yCd = [yCd ; Y(newCandidates) + Ly];
        zCd = [zCd ; Z(newCandidates) + Lz];
    end

    if xNOverlap && yPOverlap && zPOverlap
        newCandidates = find(sqrt(power(X-Lx-X(part), 2) + power(Y+Ly-Y(part), 2) + power(Z+Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) - Lx];
        yCd = [yCd ; Y(newCandidates) + Ly];
        zCd = [zCd ; Z(newCandidates) + Lz];
    end

    if xPOverlap && yNOverlap && zPOverlap
        newCandidates = find(sqrt(power(X+Lx-X(part), 2) + power(Y-Ly-Y(part), 2) + power(Z+Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) + Lx];
        yCd = [yCd ; Y(newCandidates) - Ly];
        zCd = [zCd ; Z(newCandidates) + Lz];
    end

    if xNOverlap && yNOverlap && zPOverlap
        newCandidates = find(sqrt(power(X-Lx-X(part), 2) + power(Y-Ly-Y(part), 2) + power(Z+Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) - Lx];
        yCd = [yCd ; Y(newCandidates) - Ly];
        zCd = [zCd ; Z(newCandidates) + Lz];
    end

    if xPOverlap && yPOverlap && zNOverlap
        newCandidates = find(sqrt(power(X+Lx-X(part), 2) + power(Y+Ly-Y(part), 2) + power(Z-Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) + Lx];
        yCd = [yCd ; Y(newCandidates) + Ly];
        zCd = [zCd ; Z(newCandidates) - Lz];
    end

    if xNOverlap && yPOverlap && zNOverlap
        newCandidates = find(sqrt(power(X-Lx-X(part), 2) + power(Y+Ly-Y(part), 2) + power(Z-Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) - Lx];
        yCd = [yCd ; Y(newCandidates) + Ly];
        zCd = [zCd ; Z(newCandidates) - Lz];
    end

    if xPOverlap && yNOverlap && zNOverlap
        newCandidates = find(sqrt(power(X+Lx-X(part), 2) + power(Y-Ly-Y(part), 2) + power(Z-Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) + Lx];
        yCd = [yCd ; Y(newCandidates) - Ly];
        zCd = [zCd ; Z(newCandidates) - Lz];
    end

    if xNOverlap && yNOverlap && zNOverlap
        newCandidates = find(sqrt(power(X-Lx-X(part), 2) + power(Y-Ly-Y(part), 2) + power(Z-Lz-Z(part), 2)) < 2*maxR);
        cd = [cd; newCandidates];
        xCd = [xCd ; X(newCandidates) - Lx];
        yCd = [yCd ; Y(newCandidates) - Ly];
        zCd = [zCd ; Z(newCandidates) - Lz];
    end
    
    alphaCd = alpha(cd);
    betaCd = beta(cd);
    gammaCd = gamma(cd);

end