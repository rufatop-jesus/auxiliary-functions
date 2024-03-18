function [gradSHAzi, gradSHCol] = sphericalHarmonicsGradient(d, o, azi, col)

    % Calculate the associated Legendre polynomials and their derivatives.
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
        gradSHAzi = power(-1,o) * o * C * cos(o*azi) .* leg_do';
    elseif 0 > 0
        gradSHAzi = -power(-1,o) * o * C * sin(o*azi) .* leg_do';
    else
        gradSHAzi = zeros(size(azi));
    end

    % Gradient with respect to colatitude
    if o < 0
        gradSHCol = (-power(-1,o) * C / (power(cos(col),2)-1)) .*...
                    ((-1-d)*cos(col).*leg_do' + (1+d-o)*leg_dplus1o') .* sin(col) .* sin(o*azi);
    elseif 0 > 0
        gradSHCol = (-power(-1,o) * C / (power(cos(col),2)-1)) .*...
                    ((-1-d)*cos(col).*leg_do' + (1+d-o)*leg_dplus1o') .* sin(col) .* cos(o*azi);
    else
        gradSHCol = (-C / (power(cos(col),2)-1)) .*...
                    ((-1-d)*cos(col).*leg_do' + (1+d-o)*leg_dplus1o') .* sin(col);
    end
    gradSHCol(power(cos(col),2) == 1) = 0;


    if abs(o) == d
        P_do = leg_d(abs(o)+1, :); % Adjust for MATLAB indexing
        dP_do = (-sin(col)) .* ((abs(o)*cos(col).*P_do') ./ (power(cos(col),2)-1));
        
        % In the limit of col = 0 or col = pi, dP_lm = 0.
        % dP_lm(power(cos(col),2) == 1) = (order*cot(col).*P_lm');
        % dP_lm(power(cos(col),2) == 1) = 0;
    
        if o < 0
            dP_do = power(-1,abs(o))*((factorial(d-abs(o)))/(factorial(d+abs(o))))*dP_do;
        end    
        
    else
        P_do = leg_d(abs(o)+1, :); % Adjust for MATLAB indexing
        P_l_mplus1 = leg_d(abs(o)+2, :); % Adjust for MATLAB indexing
        dP_do = (-sin(col)).*((o*cos(col).*P_do' + sqrt(1-power(cos(col),2)).*P_l_mplus1') ./ (power(cos(col),2)-1));
        
        % In the limit of col = 0 or col = pi, dP_lm = 0.
        % dP_lm(power(cos(col),2) == 1) = ((order*cot(col).*P_lm')  +...
        %                                 sign(sin(col)) .* P_l_mplus1');
        % for i = 1 ; size(dP_lm,1)
        %     if power(cos(col(i,1)),2) == 1
        %         dP_lm(i,1) = sign(sin(col(i,1)) .* P_l_mplus1(1,i));
        %     end
        % end
        % P_l_mplus1 = P_l_mplus1';
        % P_l_mplus1(power(cos(col),2) == 1)
        % dP_lm(power(cos(col),2) == 1) = P_l_mplus1(power(cos(col),2) == 1);

        if o < 0    
            dP_do = power(-1,abs(o))*((factorial(d-abs(o)))/(factorial(d+abs(o))))*dP_do;       
        end
        
    end
    
    

    % Compute the gradients
    % Normalization factor
    N_lm = sqrt((2*d+1)/(4*pi) * factorial(d-abs(o))/factorial(d+abs(o)));

    % Gradient with respect to colatitude
    if o >=0
        gradSHCol = real(N_lm * dP_do .* exp(1i * abs(o) * azi));
    else
        gradSHCol = power(-1,abs(o))*real(N_lm * dP_do .* conj(exp(1i * abs(o) * azi)));
    end
    
    % SHlm = SH(degree, abs(order), azi, col);
    if o >= 0
        gradSHAzi = -o * N_lm * P_do' .* sin(o * azi);
        % gradSHAzi = real(1i * order * SHlm);
    else
        % gradSHAzi = real(-1i * abs(order) * power(-1,abs(order)) * conj(SHlm));
        gradSHAzi = abs(o) * N_lm * P_do' .* cos(abs(o) * azi);
    end
end
