% Takes 2 solids A and B centered at
% positions OA and OB, rotated by the Euler angles alphaA, betaA, gammaA and aphaB, betaB, gammaB, and whose shape is defined according to Brechbühler,
% Gerig and Kübler (1995) by the coefficients cA and cB, and determines if
% there is contact (or overlap) between them.

% Such contact function was first developed by Garboczi and Bullard (2013).

% The uniform polar mesh in which the spherical harmonics will be
% evaluated is given by the azimuth (azi) and colaitude (col) arrays.

% Example

% % loading SH coefficients to use as example
% clear all
% addpath ("auxFunctions")
% c = readmatrix("irregular.csv","Range", 2);
% 
% % Creating a uniform polar mesh to use as example
% 
% meshRes = 100; % uniform mesh resolution    
% [uniformSphMesh,~] = spheretri(meshRes);
% 
% [azi,ele,~] = cart2sph(uniformSphMesh(:,1), uniformSphMesh(:,2), uniformSphMesh(:,3));
% 
% col = -ele + pi/2;
% for i = 1 : size(azi,1)
%     if azi(i,1) < 0
%         azi(i,1) = 2*pi + azi(i,1);
%     end
% end
% 
% OA = [0;0;0];
% OB = [1.5;0;0];
% 
% alphaA = 0;
% betaA = 0;
% gammaA = 0;
% 
% alphaB = 0;
% betaB = 0;
% gammaB = 0;
% 
% cA = c;
% cB = c;
% 
% o = Overlap(OA,alphaA,betaA,gammaA,cA, OB,alphaB,betaB,gammaB,cB, azi, col);

function o = overlap(OA,alphaA,betaA,gammaA,cA, OB,alphaB,betaB,gammaB,cB, azi, col)
    o = false;
    
    % sA = number of spherical harmonics used to describe the solid A
    sA = size(cA,2);
    
    % Calculating distances and polar coordinates of surface points of solid A (aziA,
    % eleA, rA) with respect to the center of solid B (OB)
    xBA = zeros(size(azi));
    yBA = zeros(size(azi));
    zBA = zeros(size(azi));

    degree = 1;
    order = -1;

    for i = 1 : sA  
        Y = SH(degree,order,azi,col); % real spherical harmonic function

        xBA = xBA + cA(1,i)*Y;
        yBA = yBA + cA(2,i)*Y;
        zBA = zBA + cA(3,i)*Y;
                
        % Update degree and order
        if order+1 > degree
            degree = degree+1;
            order = -degree;
        else
            order = order+1;
        end
    end

    v = [xBA,yBA,zBA];
    v = rotx(alphaA*180/pi) * roty(betaA*180/pi) * rotz(gammaA*180/pi) * v';
    
    xBA = v(1,:)' + OA(1,1) - OB(1,1);
    yBA = v(2,:)' + OA(2,1) - OB(2,1);
    zBA = v(3,:)' + OA(3,1) - OB(3,1);
    
    [aziBA,eleBA,rBA] = cart2sph(xBA,yBA,zBA);
    
    % sB = number of spherical harmonics used to describe the solid B
    sB = size(cB,2);

    % Calculating distances and polar coordinates of surface points of solid B (aziB,
    % eleB, rB) with respect to the center of solid B (OB).
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

    scatter(aziB,eleB)
    hold on
    scatter(aziBA,eleBA)
    hold off
    figure

    interprB = scatteredInterpolant(aziB,eleB,rB, 'natural','linear');
    
    % If the distance between a point of A and OB is lower than the
    % distance between the point of B and OB in the same polar coordinate,
    % then there is a collision
    if nnz(rBA < interprB(aziBA,eleBA)) > 0
        o = true;
    end

    o
    scatter3(xBA,yBA,zBA)
    hold on
    scatter3(xB,yB,zB)
    hold off
end