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
% meshRes = 100; % number of points in the uniform mesh
% PB = createPB(c,meshRes);
% [azi,ele,r] = cart2sph(PB(:,1), PB(:,2), PB(:,3)); % polar cooerdinates
% 
% % Generating a table of azimuth, elevation and radius of the
% % solid shape. The method emplys an evenly spaced grid over the sphere.
% surfaceInterpolant = scatteredInterpolant(azi,ele,r, 'linear','linear'); % Matlab interpolant function
% [aziTable,eleTable] = meshgrid((-180:1:179)*pi/180,(-90:1:89)*pi/180); % evenly spaced grid of azimuth and elevation
% rI = surfaceInterpolant(aziTable,eleTable); % interpolated radius book
% 
% OA = [0;0;0];
% OB = [0.5;0;0];
% 
% alphaA = pi/4;
% betaA = pi/6;
% gammaA = 0;
% 
% alphaB = 0;
% betaB = 0;
% gammaB = 0;
% 
% 
% o = Overlap(aziTable, eleTable, OA,alphaA,betaA,gammaA,rI, OB,alphaB,betaB,gammaB,PB);


function o = overlap(OA,alphaA,betaA,gammaA,rA, OB,alphaB,betaB,gammaB,PBB)
    o = false;

    % Definition of particle B cartesian coordinates 
    xB = PBB(:,1);
    yB = PBB(:,2);
    zB = PBB(:,3);
   
    % Rotation of particle B according to alphaB, betaB and gammaB
    v = rotx(rad2deg(alphaB)) * roty(rad2deg(betaB)) * rotz(rad2deg(gammaB)) * [xB,yB,zB]';
    
    % Translation of particle B according to OB
    xB = v(1,:) + OB(1,1) - OA(1,1);
    yB = v(2,:) + OB(2,1) - OA(2,1);
    zB = v(3,:) + OB(3,1) - OA(3,1);

    % Rotation of particle A.
    
    v = rotx(rad2deg(-alphaA)) * roty(rad2deg(-betaA)) * rotz(rad2deg(-gammaA)) * [xB;yB;zB];
 
    % Calculation of the polar coordinates of particle B relative to the
    % center of particle A.

    [aziB,eleB,rB] = cart2sph(v(1,:), v(2,:), v(3,:));

    % scatter(aziB,eleB)
    % hold on
    % scatter(aziBA,eleBA)
    % hold off
    % figure
    
    % If the distance between a point of A and OB is lower than the
    % distance between the point of B and OB in the same polar coordinate,
    % then there is a collision
    posAzi = floor(rad2deg(aziB));
    posEle = floor(rad2deg(eleB));

    for i = 1 : size(rB,2) 
        % posAzi(1,i)
        % posEle(1,i)
        if rB(1,i) < rA(posEle(1,i)+91, posAzi(1,i)+181)
            o = true;
            break
        end
    end

    % o
    % [xA,yA,zA] = sph2cart(reshape(aziTable,1,numel(aziTable)), reshape(eleTable,1,numel(eleTable)), reshape(rA,1,numel(rA)));
    % xB = v(1,:);
    % yB = v(2,:);
    % zB = v(3,:);
    % 
    % scatter3(xA,yA,zA)
    % xlim([-2 2])
    % ylim([-2 2])
    % zlim([-2 2])
    % hold on
    % scatter3(xB,yB,zB)
    % hold off
end