% Calculate the distance between particle centers and particle surfaces
% from systems with multiples particles.

% Calculate and plot frequency distributions of such distances.

function [center2center, face2face] = calculateFCCDistances(X,Y,Z,alpha,beta,gamma, rI, PB)
    nPart = size(X,1);
    center2center = c2cNearestNeighborMatrix(nPart,X,Y,Z);
    center2center = sort(center2center);

    % figure
    % plot(center2center, 0:1/(nPart-1):1)    

    face2face = f2fNearestNeighborMatrix(nPart,rI,PB,X,Y,Z,alpha,beta,gamma);
    face2face = sort(face2face);

    % figure
    % plot(face2face, 0:1/(nPart-1):1)
end

% Build the matrix with distances between particle centers
function center2center = c2cNearestNeighborMatrix(nPart,X,Y,Z)
    dMatrix = zeros(nPart);   
    for i = 1 : (nPart-1)
        for j = (i+1) : nPart
            
            % distance between particles i and j
            r = sqrt(power(X(j)-X(i), 2) + power(Y(j)-Y(i), 2) + power(Z(j)-Z(i), 2));
            dMatrix(i,j) = r;
        end
    end

    center2center = zeros(nPart,1);
    for part = 1 : nPart
        if part == 1
            center2center(part,1) = min(dMatrix(1, 2:1:nPart));
        elseif part == nPart
            center2center(part,1) = min(dMatrix(1:1:nPart-1, nPart));
        else
            center2center(part,1) = min(min(dMatrix(part, part+1:1:nPart)), min(dMatrix(1:1:part-1, part)));
        end
    end
end

function face2face = f2fNearestNeighborMatrix(nPart,rI,PB,X,Y,Z,alpha,beta,gamma)
    dMatrix = zeros(nPart);   
    for i = 1 : (nPart-1)
        for j = (i+1) : nPart
                   
            c2cvector = [X(j)-X(i); Y(j)-Y(i); Z(j)-Z(i)];

            dMatrix(i,j) = f2fDistance(c2cvector, alpha(i,1),beta(i,1),gamma(i,1),rI, alpha(j,1),beta(j,1),gamma(j,1),PB);
            if dMatrix(i,j) < 0
                dMatrix(i,j) = 0;
                % [dMatrix(i,j), i, j, X(i), Y(i), Z(i), X(j), Y(j), Z(j)]
                % plotOverlapParticles(PB,0.6567,Lx,Ly,Lz, X,Y,Z,alpha,beta,gamma, i,j)
            end
        end
    end
    face2face = zeros(nPart,1);
    for part = 1 : nPart
        if part == 1
            face2face(part,1) = min(dMatrix(1, 2:1:nPart));
        elseif part == nPart
            face2face(part,1) = min(dMatrix(1:1:nPart-1, nPart));
        else
            face2face(part,1) = min(min(dMatrix(part, part+1:1:nPart)), min(dMatrix(1:1:part-1, part)));
        end
    end
end

function d = f2fDistance(c2cvector, alphaA,betaA,gammaA,rA, alphaB,betaB,gammaB,PBB)
    dc2c = sqrt(power(c2cvector(1,1),2) + power(c2cvector(2,1),2) + power(c2cvector(3,1),2));
    d = dc2c;
    
    % Definition of particle B cartesian coordinates 
    xB = PBB(:,1);
    yB = PBB(:,2);
    zB = PBB(:,3);
   
    % Rotation of particle B according to alphaB, betaB and gammaB
    v = rotx(rad2deg(alphaB)) * roty(rad2deg(betaB)) * rotz(rad2deg(gammaB)) * [xB,yB,zB]';
    
    % Translation of particle B according to the c2cvector
    xB = v(1,:) + c2cvector(1,1);
    yB = v(2,:) + c2cvector(2,1);
    zB = v(3,:) + c2cvector(3,1);

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
        
        if d > rB(1,i) - rA(posEle(1,i)+91, posAzi(1,i)+181)
            d = rB(1,i) - rA(posEle(1,i)+91, posAzi(1,i)+181);
        end
    end
end

function S = plotOverlapParticles(PB,maxR,Lx,Ly,Lz, x,y,z,alpha,beta,gamma, A,B)

    % Rotating and translating the shape A
    v = rotx(rad2deg(alpha(A))) * roty(rad2deg(beta(A))) * rotz(rad2deg(gamma(A))) * PB';
    X = v(1,:)' + x(A);
    Y = v(2,:)' + y(A);
    Z = v(3,:)' + z(A);
    shpA = alphaShape(X,Y,Z,maxR);

    plot(shpA,'FaceColor',[0.6350 0.0780 0.1840],'EdgeAlpha',0)
    hold on

    % Rotating and translating the shape B
    v = rotx(rad2deg(alpha(B))) * roty(rad2deg(beta(B))) * rotz(rad2deg(gamma(B))) * PB';
    X = v(1,:)' + x(B);
    Y = v(2,:)' + y(B);
    Z = v(3,:)' + z(B);
    shpB = alphaShape(X,Y,Z,maxR);

    plot(shpB,'FaceColor',[0.6350 0.0780 0.1840],'EdgeAlpha',0)
    hold on

    lighting gouraud
    light('Position',[2 -4 2],'Style','local')
    material dull
    xlim([0 Lx])
    ylim([0 Ly])
    zlim([0 Lz])
    hold off
    figure
end