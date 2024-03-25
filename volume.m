% Calculate the volume of a particle defined by the Brechbühler, Gerig and Kübler (1995)
% method using spherical harmonics.

% clear all
% addpath("auxFunctions")
% c = readmatrix("irregular.csv","Range", 2);
% meshRes = 10; % number of points on the solid surface
% [uniformSphMesh,~] = spheretri(meshRes); % creating a uniform spherical mesh
% uniformSphMesh = transpose(rotx(90*rand)*roty(90*rand)*rotz(90*rand)*uniformSphMesh'); % random rotation to avoid points on the North and South poles.
% [azi,ele,~] = cart2sph(uniformSphMesh(:,1), uniformSphMesh(:,2), uniformSphMesh(:,3)); % transforming from caetesian to spherical coordinates
% col = -ele + pi/2;
% for i = 1 : size(azi,1)
%     if azi(i,1) < 0
%         azi(i,1) = 2*pi + azi(i,1);
%     end
% end
% v = Volume(c,azi,col);

function v = volume(c)
    s = size(c,2);

    % UNIFORM MESH
    meshRes = 1000; % number of points on the solid surface
    [uniformSphMesh,~] = spheretri(meshRes); % creating a uniform spherical mesh
    uniformSphMesh = transpose(rotx(90*rand)*roty(90*rand)*rotz(90*rand)*uniformSphMesh'); % random rotation to avoid points on the North and South poles.
    [azi,ele,~] = cart2sph(uniformSphMesh(:,1), uniformSphMesh(:,2), uniformSphMesh(:,3)); % transforming from caetesian to spherical coordinates
    col = -ele + pi/2;
    for i = 1 : size(azi,1)
        if azi(i,1) < 0
            azi(i,1) = 2*pi + azi(i,1);
        end
    end

    % Declaring gradX, gradY, and gradZ vectors   
    gradX = zeros(size(azi,1),2);
    gradY = zeros(size(azi,1),2);
    gradZ = zeros(size(azi,1),2);

    x = zeros(size(azi,1),1);
    y = zeros(size(azi,1),1);
    z = zeros(size(azi,1),1);

    d = 1; % order
    o = -1; % degree

    for i = 1 : s  
        [gradSHAzi, gradSHCol] = sphericalHarmonicsGradient(d, o, azi, col); % the real spherical harmonic gradient function
        Y = SH(d, o, azi, col);
        
        % Feeding gradX, gradY and gradZ values with SHGradient values according to the
        % respective coefficients in matrix c
        gradX = gradX + c(1,i)*[gradSHAzi, gradSHCol];
        gradY = gradY + c(2,i)*[gradSHAzi, gradSHCol];
        gradZ = gradZ + c(3,i)*[gradSHAzi, gradSHCol];

        % Feeding x, y and z values with SH values according to the
        % respective coefficients in matrix c
        x = x + c(1,i)*Y;
        y = y + c(2,i)*Y;
        z = z + c(3,i)*Y;
                
        % Update degree and order
        if o+1 > d
            d = d+1;
            o = -d;
        else
            o = o+1;
        end
    end
    
    % Calculation of the radius
    [~,~,r] = cart2sph(x,y,z);

    dS = vecnorm(cross([gradX(:,1),gradY(:,1),gradZ(:,1)],[gradX(:,2),gradY(:,2),gradZ(:,2)]), 2, 2);
    dS = dS ./ sin(col);
    dO = 4*pi / size(azi,1); % dO is the average solid angle of the surface uniform mesh
    dS = dS * dO;

    dV = r .* dS / 3;

    v = sum(dV);
end