% Calculate the volume of a particle defined by the Brechbühler, Gerig and Kübler (1995)
% method using spherical harmonics.

% clear all
% addpath("auxFunctions")
% c = readmatrix("sphere.csv","Range", 2);
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
% [xmax, xmin, ymax, ymin, zmax, zmin] = Box(c);

function [xFeret, yFeret, zFeret] = box(c)
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

    x = zeros(size(azi,1),1);
    y = zeros(size(azi,1),1);
    z = zeros(size(azi,1),1);

    d = 1; % order
    o = -1; % degree

    for i = 1 : s  
        Y = SH(d, o, azi, col);

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
    
    xFeret = max(x) - min(x);
    yFeret = max(y) - min(y);
    zFeret = max(z) - min(z);


    % xmax = max(x);
    % xmin = min(x);
    % ymax = max(y);
    % ymin = min(y);
    % zmax = max(z);
    % zmin = min(z);
end