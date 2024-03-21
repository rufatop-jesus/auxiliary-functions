clear all
addpath ("auxFunctions")

% Takes 2 solids A and B centered at
% positions OA and OB, rotated by the Euler angles alphaA, betaA, gammaA and aphaB, betaB, gammaB, and whose shape is defined according to Brechbühler,
% Gerig and Kübler (1995) by the coefficients cA and cB, and determines if
% there is contact (or overlap) between them.

% Such contact function was first developed by Garboczi and Bullard (2013).

% The uniform polar mesh in which the spherical harmonics will be
% evaluated is given by the azimuth (azi) and colaitude (col) arrays.

% loading SH coefficients to use as example
c = readmatrix("sphere.csv","Range", 2);

% Creating a uniform polar mesh to use as example

meshRes = 10; % uniform mesh resolution    
[uniformSphMesh,~] = spheretri(meshRes);

[azi,ele,~] = cart2sph(uniformSphMesh(:,1), uniformSphMesh(:,2), uniformSphMesh(:,3));

col = -ele + pi/2;
for i = 1 : size(azi,1)
    if azi(i,1) < 0
        azi(i,1) = 2*pi + azi(i,1);
    end
end

OA = [0;0;0];
OB = [1;0;0];

alphaA = 0;
betaA = 0;
gammaA = 0;

alphaB = 0;
betaB = 0;
gammaB = 0;

cA = c;
cB = c;

o = Overlap(OA,alphaA,betaA,gammaA,cA, OB,alphaB,betaB,gammaB,cB, azi, col);

function o = Overlap(OA,alphaA,betaA,gammaA,cA, OB,alphaB,betaB,gammaB,cB, azi, col)
    o = false;
    
    % sA = number of spherical harmonics used to describe the solid A
    sA = size(cA,2);
    
    % Calculating distances and polar coordinates of surface points of solid A (aziA,
    % eleA, rA) with respect to the origin [0;0;0]
    xA = zeros(size(azi));
    yA = zeros(size(azi));
    zA = zeros(size(azi));

    degree = 1;
    order = -1;

    for i = 1 : sA  
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
    
    % [aziA,eleA,rA] = cart2sph(xA,yA,zA);
    
    % sB = number of spherical harmonics used to describe the solid B
    sB = size(cB,2);

    % Calculating distances and polar coordinates of surface points of solid B (aziB,
    % eleB, rB) with respect to the center of B (OB).
    xB = zeros(size(azi));
    yB = zeros(size(azi));
    zB = zeros(size(azi));

    degree = 1;
    order = -1;

    for i = 1 : sB  
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
    
    xB = v(1,:)';
    yB = v(2,:)';
    zB = v(3,:)';
    
    [aziB,eleB,rB] = cart2sph(xB,yB,zB);

    % Calculating the box around particles A and B
    % boxA = [min(xA),max(xA);...
    %         min(xA),max(xA);...
    %         min(zA),max(zA)];
    % 
    % boxB = [min(xB),max(xB);...
    %         min(xB),max(xB);...
    %         min(zB),max(zB)];

    % Calculating distances and polar coordinates between the surface points of A to the center of
    % B (OB)
    xBA = xA - OB(1,1);
    yBA = yA - OB(2,1);
    zBA = zA - OB(3,1);

    [aziBA,eleBA,rBA] = cart2sph(xBA,yBA,zBA);

    scatter(aziB,eleB)
    hold on
    scatter(aziBA,eleBA)
    hold off
    figure

    interprB = griddata(aziB,eleB,rB, aziBA,eleBA);
    
    % If the distance between a point of A and OB is lower than the
    % distance between the point of B and OB in the same polar coordinate,
    % then there is a collision
    if nnz((rBA - interprB) > 0) > 0
        o = true;
    end
    
    o
    scatter3(xA,yA,zA)
    hold on
    scatter3(xB + OB(1,:),yB + OB(2,:),zB + OB(3,:))
    hold off
end