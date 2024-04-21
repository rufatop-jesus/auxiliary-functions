% Identify which walls are overlapped, if any.
function [xPOverlap,xNOverlap, yPOverlap,yNOverlap, zPOverlap,zNOverlap] = identifyWall(maxR, Lx,Ly,Lz, X,Y,Z)
    xPOverlap = false;
    xNOverlap = false;

    yPOverlap = false;
    yNOverlap = false;

    zPOverlap = false;
    zNOverlap = false;

    if X < 2*maxR         
        xNOverlap = true;
    
    elseif X > Lx - 2*maxR
        xPOverlap = true;
    
    end

    if Y < 2*maxR         
        yNOverlap = true;
    
    elseif Y > Ly - 2*maxR
        yPOverlap = true;
    
    end

    if Z < 2*maxR         
        zNOverlap = true;
    
    elseif Z > Lz - 2*maxR
        zPOverlap = true;
    
    end
end