% Repository of the functions used to calculate shape factors

% X is the matrix containing cartesian coordinates of the 2D projection
function [figPerimeter, figArea, feretMax, feretMin] = SFScatter(X, polarRes, type)
    
    % bPC = boundaryPolarCoordinates
    bPC = pBoundary(X, polarRes,type);
        
    % Calculation of the area
    [x, y] = pol2cart(bPC(:,1), bPC(:,2));

    figArea = 0;

    for p = 1 : size(x,1)
        if p == 1
            dx = x(size(x,1),1) - x(2,1);
            dA = y(p,1)*dx;
        elseif p == size(x,1)
            dx = x(size(x,1)-1,1) - x(1,1);
            dA = y(p,1)*dx;
        else
            dx = x(p-1, 1) - x(p+1, 1);
            dA = y(p,1)*dx;
        end
        figArea = figArea + dA;
    end
    figArea = figArea/2;

    % Numerically integrating the radius data to calculate perimeter.
    % Source: https://web.ma.utexas.edu/users/m408s/m408d/CurrentWeb/LM10-4-4.php

    figPerimeter = trapz(bPC(:,1),bPC(:,4));
    
    % polarplot(bPC(:,1), bPC(:,2))
    

    % Calculating the Feret diameters. WORKS BEST FOR STAR-LIKE PROJECTIONS
    feretDiameter = zeros(180, 1);
    for k = 0 : 179 
        feretDiameter(k+1,1) = interp1(bPC(:,1), bPC(:,2), k*pi/180) + interp1(bPC(:,1), bPC(:,2), (k-180)*pi/180);
    end
    feretMax = maxk(feretDiameter,1);
    feretMin = mink(feretDiameter,1);
    
end

% input: cartesian coordinates of a figure, defined as X
% output: polar coordinates of the figure's boundary

function sbPC = pBoundary(X, polarRes, type)
    
    if type == "custom"
    
        % projection center of gravity
        xcg = mean(X(:,1)); 
        ycg = mean(X(:,2));
    
        s = size(X,1);
    
        if rem(polarRes, 2) ~= 0
            polarRes = polarRes - 1;
        end
        dTheta = 2*pi/polarRes; % angle differential
        
        % Identifying the boundaries 
        sbPC = zeros(polarRes,4); % sorted boundary polar coordinates. (angle, r, r', sqrt(r + r')) 
            
        % Building the angles array
        for k = 1 : polarRes
            sbPC(k,1) = -pi + (k-1)*dTheta;
        end
        
        % Calculating polar from cartesian coordinates
        for k = 1 : s
            radTest = sqrt(power(X(k,1)-xcg,2)+power(X(k,2)-ycg,2));
    
            if(X(k,1)-xcg > 0)
    
                angle = atan( (X(k,2)-ycg) / (X(k,1)-xcg) );
            
            elseif(X(k,1)-xcg < 0)
    
                angle = atan( (X(k,2)-ycg) / (X(k,1)-xcg) ) + pi*sign(X(k,2)-ycg);
    
            else % if(XP(k,1)-xcg == 0)
    
                if(X(k,2)-ycg > 0)
    
                    angle = pi/2;
    
                elseif(X(k,2)-ycg < 0)
    
                    angle = -pi/2;
    
                else
                    angle = 0;
                end
            end

            if (angle > -pi)
                position = ceil((angle+pi)/dTheta);
            elseif(angle == -pi)
                position = polarRes;
            end
            if sbPC (position,2) < radTest
                sbPC(position,2) = radTest;
            end
    
           sbPC(:,3) = gradient(sbPC(:,2), dTheta); % radius derivative with respect to the angle
           sbPC(:,4) = sqrt(power(sbPC(:,2),2) + power(sbPC(:,3),2)); % integrand for the calculus of the perimeter      
        end  
    elseif type == "built-in"
        
        % Identifying the boundaries        
        boundaryPos = boundary(X,1);
        boundaryCartCoordinates = X(boundaryPos,:);

        % Max radius list. The angle in degrees is a multiple of the row number. (angle, r, r^2)
        [angle, radius] = cart2pol(boundaryCartCoordinates(:,1), boundaryCartCoordinates(:,2));
        
        sbPC = zeros(size(boundaryCartCoordinates,1),4);
        sbPC(:,1) = angle; 
        sbPC(:,2) = radius;
        sbPC(:,3) = gradient(sbPC(:,2), sbPC(:,1)); % radius derivative with respect to the angle
        sbPC(:,4) = sqrt(power(sbPC(:,2),2) + power(sbPC(:,3),2)); % integrand for the calculus of the perimeter
        
        % Eliminating repeated points. ubPC is the unique version of
        % boundaryPolarCoordinates
        [~,uniquebPC] = unique(sbPC(:,1));
        ubPC = [sbPC(uniquebPC,1), sbPC(uniquebPC,2),...
               sbPC(uniquebPC,3),sbPC(uniquebPC,4)];

        % sbPC is the sorted version of boundaryPolarCoordinates
        sbPC = zeros(size(ubPC));       
        [sbPC(:,1), I] = sort(ubPC(:,1));

        sbPC(:,2) = ubPC(I,2);
        sbPC(:,3) = ubPC(I,3);
        sbPC(:,4) = ubPC(I,4);
    end
end