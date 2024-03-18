function [gradSHAzi, gradSHCol] = sphericalHarmonicsGradient(d, o, azi, col)

    % Calculate the spherical harmonics derivatives with respect to azimuth
    % and colatitude.

    % degree -> d, order -> o
    
    % Normalization factor
    a = (2*d+1)*factorial(d-abs(o));
    b = 4*pi*factorial(d+abs(o));
    C = sqrt(a/b);

    % Legendre polynomials
    leg_d = legendre(d, cos(col));
    if d > 0
        leg_do = leg_d(abs(o)+1, :);
    else
        leg_do = leg_d(1, :);
    end

    leg_dplus1 = legendre(d+1, cos(col));
    if d > 0
        leg_dplus1o = leg_dplus1(abs(o)+1, :);
    else
        leg_dplus1o = leg_dplus1(1, :);
    end
    
    % Gradient with respect to azimuth
    if o < 0
        gradSHAzi = sqrt(2) * power(-1,o) * abs(o) * C * cos(abs(o)*azi) .* leg_do';
    elseif o > 0
        gradSHAzi = -sqrt(2) * power(-1,o) * o * C * sin(o*azi) .* leg_do';
    else
        gradSHAzi = zeros(size(azi));
    end

    % Gradient with respect to colatitude
    if o < 0
        gradSHCol = -sqrt(2) * power(-1,o) * C * ((-1-d)*cos(col).*leg_do' + (1+d-abs(o))*leg_dplus1o') .*....
                    sin(col) .* sin(abs(o)*azi) ./ (power(cos(col),2)-1);
    elseif o > 0
        gradSHCol = -sqrt(2) * power(-1,o) * C * ((-1-d)*cos(col).*leg_do' + (1+d-o)*leg_dplus1o') .*...
                    sin(col) .* cos(o*azi) ./ (power(cos(col),2)-1);
    else
        gradSHCol = -C * ((-1-d)*cos(col).*leg_do' + (1+d-o)*leg_dplus1o') .*...
                     sin(col) ./ (power(cos(col),2)-1);
    end
    % gradSHCol(power(cos(col),2) == 1) = 0;
    
    % gradSHAzi
    % gradSHCol
end
