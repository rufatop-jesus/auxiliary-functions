% input: cartesian coordinates of a figure, defined as X
% output: polar coordinates of the figure's boundary

function boundaryPolarCoordinates = customBoundary(X, polarRes)

    % projection center of gravity
    xcg = mean(X(:,1)); 
    ycg = mean(X(:,2));

    s = size(X,1);

    if rem(polarRes, 2) ~= 0
        polarRes = polarRes - 1;
    end
    dTheta = 2*pi/polarRes; % angle differential
    
    % Identifying the boundaries 
    boundaryPolarCoordinates = zeros(polarRes,4); % max radius list. The angle in degrees is a multiple of the row number. (angle, r, r^2) 
        
    % Building the angles array
    for k = 1 : polarRes
        boundaryPolarCoordinates(k,1) = k*dTheta;
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
        elseif(angle == 0)
            position = polarRes;
        end
        if boundaryPolarCoordinates (position,2) < radTest
            boundaryPolarCoordinates(position,2) = radTest;
        end
    end

    boundaryPolarCoordinates(:,3) = gradient(boundaryPolarCoordinates(:,2), dTheta); % radius derivative with respect to the angle
    boundaryPolarCoordinates(:,4) = sqrt(power(boundaryPolarCoordinates(:,2),2) + power(boundaryPolarCoordinates(:,3),2)); % integrand for the calculus of the perimeter
end