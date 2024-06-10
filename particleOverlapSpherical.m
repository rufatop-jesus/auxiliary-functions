% Verify if there is particle overlap between the particle "part" and all
% the candidates in "cd" with position x,y,z - cd and rotation alpha,beta,gamma cd.
function o = particleOverlapSpherical(collisionCandidates,cellSize,part,partDiameter,L,X,Y,Z)
    o = false;

    for i = collisionCandidates

        if i == part
            continue
        end
        
        % dx = min([abs(X(i)-X(part)), abs(X(i)-X(part)+L), abs(X(i)-X(part)-L)]);
        % dy = min([abs(Y(i)-Y(part)), abs(Y(i)-Y(part)+L), abs(Y(i)-Y(part)-L)]);
        % dz = min([abs(Z(i)-Z(part)), abs(Z(i)-Z(part)+L), abs(Z(i)-Z(part)-L)]);

        if abs(X(part)-X(i)) > partDiameter + cellSize
            if X(i) < X(part)
                dx = X(i) - X(part) + L;
            else
                dx = X(i) - X(part) - L;
            end
        else
            dx = X(i) - X(part);
        end

        if abs(Y(part)-Y(i)) > partDiameter + cellSize
            if Y(i) < Y(part)
                dy = Y(i) - Y(part) + L;
            else
                dy = Y(i) - Y(part) - L;
            end
        else
            dy = Y(i) - Y(part);
        end

        if abs(Z(part)-Z(i)) > partDiameter + cellSize
            if Z(i) < Z(part)
                dz = Z(i) - Z(part) + L;
            else
                dz = Z(i) - Z(part) - L;
            end
        else
            dz = Z(i) - Z(part);
        end

        dist = sqrt(power(dx, 2) + power(dy, 2) + power(dz, 2));

        % distance between particles i and part in the domain and in every
        % 26 periodic boundaries
        % dist = [sqrt(power(X(part)-X(i), 2) + power(Y(part)-Y(i), 2) + power(Z(part)-Z(i), 2));...
        % 
        %         sqrt(power(X(part)-X(i)+L, 2) + power(Y(part)-Y(i), 2) + power(Z(part)-Z(i), 2));...
        %         sqrt(power(X(part)-X(i)-L, 2) + power(Y(part)-Y(i), 2) + power(Z(part)-Z(i), 2));...
        % 
        %         sqrt(power(X(part)-X(i), 2) + power(Y(part)-Y(i)+L, 2) + power(Z(part)-Z(i), 2));...
        %         sqrt(power(X(part)-X(i), 2) + power(Y(part)-Y(i)-L, 2) + power(Z(part)-Z(i), 2));...
        % 
        %         sqrt(power(X(part)-X(i), 2) + power(Y(part)-Y(i), 2) + power(Z(part)-Z(i)+L, 2));...
        %         sqrt(power(X(part)-X(i), 2) + power(Y(part)-Y(i), 2) + power(Z(part)-Z(i)-L, 2));...
        % 
        %         sqrt(power(X(part)-X(i)+L, 2) + power(Y(part)-Y(i)+L, 2) + power(Z(part)-Z(i), 2));...
        %         sqrt(power(X(part)-X(i)+L, 2) + power(Y(part)-Y(i)-L, 2) + power(Z(part)-Z(i), 2));...
        %         sqrt(power(X(part)-X(i)-L, 2) + power(Y(part)-Y(i)+L, 2) + power(Z(part)-Z(i), 2));...
        %         sqrt(power(X(part)-X(i)-L, 2) + power(Y(part)-Y(i)-L, 2) + power(Z(part)-Z(i), 2));...
        % 
        %         sqrt(power(X(part)-X(i)+L, 2) + power(Y(part)-Y(i), 2) + power(Z(part)-Z(i)+L, 2));...
        %         sqrt(power(X(part)-X(i)+L, 2) + power(Y(part)-Y(i), 2) + power(Z(part)-Z(i)-L, 2));...
        %         sqrt(power(X(part)-X(i)-L, 2) + power(Y(part)-Y(i), 2) + power(Z(part)-Z(i)+L, 2));...
        %         sqrt(power(X(part)-X(i)-L, 2) + power(Y(part)-Y(i), 2) + power(Z(part)-Z(i)-L, 2));...
        % 
        %         sqrt(power(X(part)-X(i), 2) + power(Y(part)-Y(i)+L, 2) + power(Z(part)-Z(i)+L, 2));...
        %         sqrt(power(X(part)-X(i), 2) + power(Y(part)-Y(i)+L, 2) + power(Z(part)-Z(i)-L, 2));...
        %         sqrt(power(X(part)-X(i), 2) + power(Y(part)-Y(i)-L, 2) + power(Z(part)-Z(i)+L, 2));...
        %         sqrt(power(X(part)-X(i), 2) + power(Y(part)-Y(i)-L, 2) + power(Z(part)-Z(i)-L, 2));...
        % 
        %         sqrt(power(X(part)-X(i)+L, 2) + power(Y(part)-Y(i)+L, 2) + power(Z(part)-Z(i)+L, 2));...
        %         sqrt(power(X(part)-X(i)+L, 2) + power(Y(part)-Y(i)-L, 2) + power(Z(part)-Z(i)+L, 2));...
        %         sqrt(power(X(part)-X(i)-L, 2) + power(Y(part)-Y(i)+L, 2) + power(Z(part)-Z(i)+L, 2));...
        %         sqrt(power(X(part)-X(i)-L, 2) + power(Y(part)-Y(i)-L, 2) + power(Z(part)-Z(i)+L, 2));...
        % 
        %         sqrt(power(X(part)-X(i)+L, 2) + power(Y(part)-Y(i)+L, 2) + power(Z(part)-Z(i)-L, 2));...
        %         sqrt(power(X(part)-X(i)+L, 2) + power(Y(part)-Y(i)-L, 2) + power(Z(part)-Z(i)-L, 2));...
        %         sqrt(power(X(part)-X(i)-L, 2) + power(Y(part)-Y(i)+L, 2) + power(Z(part)-Z(i)-L, 2));...
        %         sqrt(power(X(part)-X(i)-L, 2) + power(Y(part)-Y(i)-L, 2) + power(Z(part)-Z(i)-L, 2))]; 

        % dist = min(dist,[],"all");

        if dist < partDiameter
            o = true;
            break
        end
    end
end