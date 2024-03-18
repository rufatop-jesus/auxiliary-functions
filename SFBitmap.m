% Calculate perimeter and Feret diameters from a figure defined by
% the list of coordinates of the points in the boundary, X.

% It is used when reading the bitmap data

function [perimeter, feretMax, feretMin] = SFBitmap(X, polarRes)
    
    % projection center of gravity
    xcg = mean(X(:,1)); 
    ycg = mean(X(:,2));

    s = size(X,1);  

    if rem(polarRes, 2) == 0
    else
        polarRes = polarRes - 1;
    end
    dTheta = 2*pi/polarRes; % angle differential
    
    % Identifying the boundaries 
    radius = zeros(polarRes,4); % max radius list. The angle in degrees is a multiple of the row number. (angle, r, r^2) 
    feretDiameter = zeros(polarRes/2, 1);
    
    % Building the angles array
    for k = 1 : polarRes
        radius(k,1) = k*dTheta;
    end
    
    % Calculating polar from cartesian coordinates
    for k = 1 : s
        radTest = sqrt(power(X(k,1)-xcg,2)+power(X(k,2)-ycg,2));

        if(X(k,1)-xcg > 0 && X(k,2)-ycg >= 0)

            angle = atan( (X(k,2)-ycg) / (X(k,1)-xcg) );
        
        elseif(X(k,1)-xcg < 0)

            angle = atan( (X(k,2)-ycg) / (X(k,1)-xcg) ) + pi;

        elseif(X(k,1)-xcg > 0 && X(k,2)-ycg < 0)

            angle = atan( (X(k,2)-ycg) / (X(k,1)-xcg) ) + 2*pi;

        else % if(XP(k,1)-xcg == 0)

            if(X(k,2)-ycg > 0)

                angle = pi/2;

            elseif(X(k,2)-ycg < 0)

                angle = 3*pi/2;

            else
                angle = 0;
            end
        end
        if (angle ~= 0)
            position = ceil(angle/dTheta);
        else
            position = 1;
        end
        if radius (position,2) < radTest
            radius(position,2) = radTest;
        end
    end

    radius(:,3) = gradient(radius(:,2), dTheta); % radius derivative with respect to the angle
    radius(:,4) = sqrt(power(radius(:,2),2) + power(radius(:,3),2)); % integrand for the calculus of the perimeter
    
    % Numerically integrating the radius data to calculate perimeter.
    % Source: https://web.ma.utexas.edu/users/m408s/m408d/CurrentWeb/LM10-4-4.php

    perimeter = trapz(dTheta,radius(:,4));
    
    % polarplot(radius(:,1), sqrt(radius(:,2)))
    % figure
    % plot(radius(:,1), sqrt(radius(:,2)))

    % Calculating the Faret diameters
    for k = 1 : polarRes/2
        feretDiameter(k,1) = radius(k,2) + radius(k+polarRes/2,2);
    end
    feretMax = maxk(feretDiameter,1);
    feretMin = mink(feretDiameter,1);
end