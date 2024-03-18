% Calculate the real spherical harmonics Ydo of some degree (d) and order (o) on a grid of spherical angles as
% defined in https://www.mathworks.com/help/matlab/ref/legendre.html and 
% https://en.wikipedia.org/wiki/Spherical_harmonics in the real form
% convention section


% TESTING VALUES, keep commented while using the function
% colatitude counts from the z-positive axis from 0 to pi. Azimuth conts from the x
% positive axis from 0 to 2*pi
% 
% azi = pi/0.7;
% col = pi/9;
% degree = 1;
% order = -1;
% Ylm = rSH(degree, order, azi, col)
% 
% if degree == 0
%     YlmAnalytic = real((1/2)*sqrt(1/pi)) % l=0, m=0
% 
% elseif degree == 1
%     if order == -1
%         YlmAnalytic = real((1/2)*sqrt(3/(2*pi))*sin(col)*exp(-1i*azi)) % l=1, m=-1
% 
%     elseif order == 0
%         YlmAnalytic = real((1/2)*sqrt(3/pi)*cos(col)) % l=1, m=0
% 
%     elseif order == 1
%         YlmAnalytic = real((-1/2)*sqrt(3/(2*pi))*sin(col)*exp(1i*azi)) % l=1, m=1
%     end
% 
% end


function Ylm = SH(d, o, azi, col)
    Plm = legendre(d,cos(col));
    if d > 0
        Plm = reshape(Plm(abs(o)+1,:,:),size(azi));
        a = (2*d+1)*factorial(d-abs(o));
        b = 4*pi*factorial(d+abs(o));
        C = sqrt(a/b);
        Ylm = C .* Plm .* exp(1i*abs(o)*azi);        
        if o < 0
            Ylm = sqrt(2)*power(-1,o)*imag(Ylm);
        elseif o > 0
            Ylm = sqrt(2)*power(-1,o)*real(Ylm);           
        end
    else
        Plm = reshape(Plm(1,:,:),size(azi));
        a = 1;
        b = 4*pi;
        C = sqrt(a/b);
        Ylm = C .* Plm;
    end
end