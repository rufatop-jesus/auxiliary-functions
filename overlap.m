clear all
addpath ("auxFunctions")

% Takes 2 solids A and B centered at
% positions OA and OB, rotated by the Euler angles alphaA, betaA, gammaA and aphaB, betaB, gammaB, and whose shape is defined according to Brechbühler,
% Gerig and Kübler (1995) by the coefficients cA and cB, and determines if
% there is contact (or overlap) between them.

% Such contact function was first developed by Garboczi and Bullard (2013).

% The uniform polar mesh in which the spherical harmonics will be
% evaluated is given by the azimuth (azi) and colaitude (col) arrays.

% Creating a uniform polar mesh to use as example
    
[uniformSphMesh,~] = spheretri(meshRes);

[azi,ele,~] = cart2sph(uniformSphMesh(:,1), uniformSphMesh(:,2), uniformSphMesh(:,3));

col = -ele + pi/2;
for i = 1 : size(azi,1)
    if azi(i,1) < 0
        azi(i,1) = 2*pi + azi(i,1);
    end
end

function o = Overlap(OA,alphaA,betaA,gammaA,cA, OB,alphaB,betaB,gammaB,cB, azi, col)
    
    % s = number of spherical harmonics used to describe the solid
    s = size(c,2);
    
    % Calculating polar coordinates of surface points of solid A (aziA,
    % eleA, rA)
    xA = zeros(size(azi));
    yA = zeros(size(azi));
    zA = zeros(size(azi));

    degree = 1;
    order = -1;

    for i = 1 : s  
        Y = SH(degree,order,azi,col); % real spherical harmonic function

        xA = xA + cA(1,i)*Y;
        yA = yA + cA(2,i)*Y;
        zA = zA + cA(3,i)*Y;
                
        % Update degree and order
        if order+1 > degree
            degree = degree+1;
            order = -degree;
        else
            order = order+1;
        end
    end

    v = [xA,yA,zA];
    v = rotx(alphaA*180/pi) * roty(betaA*180/pi) * rotz(gammaA*180/pi) * v';
    
    xA = v(1,:)' + OA(1,1);
    yA = v(2,:)' + OA(2,1);
    zA = v(3,:)' + OA(3,1);
    
    [aziA,eleA,rA] = cart2sph(xA,yA,zA);

    % Calculating polar coordinates of surface points of solid B (aziB,
    % eleB, rB)
    xB = zeros(size(azi));
    yB = zeros(size(azi));
    zB = zeros(size(azi));

    degree = 1;
    order = -1;

    for i = 1 : s  
        Y = SH(degree,order,azi,col); % real spherical harmonic function

        xB = xB + cB(1,i)*Y;
        yB = yB + cB(2,i)*Y;
        zB = zB + cB(3,i)*Y;
                
        % Update degree and order
        if order+1 > degree
            degree = degree+1;
            order = -degree;
        else
            order = order+1;
        end
    end

    v = [xB,yB,zB];
    v = rotx(alphaB*180/pi) * roty(betaB*180/pi) * rotz(gammaB*180/pi) * v';
    
    xB = v(1,:)' + OB(1,1);
    yB = v(2,:)' + OB(2,1);
    zB = v(3,:)' + OB(3,1);
    
    [aziB,eleB,rB] = cart2sph(xB,yB,zB);

    % Calculating the box around particles A and B
    boxA = [min(xA),max(xA);...
            min(xA),max(xA);...
            min(zA),max(zA)];

    boxB = [min(xB),max(xB);...
            min(xB),max(xB);...
            min(zB),max(zB)];

    
end