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
                position = 1;
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