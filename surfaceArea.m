clear all

c = readmatrix("sphere.csv","Range", 2);

% CARTESIAN MESH
% meshRes = 10; % number of points between the North and Spouth poles
% [azi, col] = meshgrid(0: pi/meshRes: pi*(2-1/meshRes), pi/meshRes: pi/meshRes: pi*(1-1/meshRes));
% azi = reshape(azi,numel(azi),1);
% col = reshape(col,numel(col),1);

% UNIFORM MESH
meshRes = 10000; % number of points on the solid surface
[uniformSphMesh,~] = spheretri(meshRes); % creating a uniform spherical mesh
uniformSphMesh = transpose(rotx(90*rand)*roty(90*rand)*rotz(90*rand)*uniformSphMesh'); % random rotation to avoid points on the North and South poles.
[azi,ele,~] = cart2sph(uniformSphMesh(:,1), uniformSphMesh(:,2), uniformSphMesh(:,3)); % transforming from caetesian to spherical coordinates
col = -ele + pi/2;
for i = 1 : size(azi,1)
    if azi(i,1) < 0
        azi(i,1) = 2*pi + azi(i,1);
    end
end

S = SA(c,azi,col,meshRes);

function S = SA(c,azi,col,meshRes)
    s = size(c,2);

    % Declaring gradX, gradY, and gradZ vectors   
    gradX = zeros(size(azi,1),2);
    gradY = zeros(size(azi,1),2);
    gradZ = zeros(size(azi,1),2);

    d = 1; % order
    o = -1; % degree

    for i = 1 : s  
        [gradSHAzi, gradSHCol] = sphericalHarmonicsGradient(d, o, azi, col); % the real spherical harmonic gradient function
        
        % Feeding gradX, gradY and gradZ values with SHGradient values according to the
        % respective coefficients in matrix c
        gradX = gradX + c(1,i)*[gradSHAzi, gradSHCol];
        % mean(abs(gradX(:,1)))
        % mean(abs(gradX(:,2)))
        gradY = gradY + c(2,i)*[gradSHAzi, gradSHCol];
        gradZ = gradZ + c(3,i)*[gradSHAzi, gradSHCol];
                
        % Update degree and order
        if o+1 > d
            d = d+1;
            o = -d;
        else
            o = o+1;
        end
    end
    % scatter3(azi,col,gradX(:,1)) % Correct
    % figure
    % scatter3(azi,col,gradY(:,1)) % Correct
    % figure
    % scatter3(azi,col,gradZ(:,1)) % Correct
    % figure
    % scatter3(azi,col,gradX(:,2)) % Correct
    % figure
    % scatter3(azi,col,gradY(:,2)) % Correct
    % figure
    % scatter3(azi,col,gradZ(:,2)) % Correct
    % figure

    % deivative of r(x,y,z) concerning azimuth
    mean(vecnorm([gradX(:,1),gradY(:,1),gradZ(:,1)], 2, 2));
    
    % deivative of r(x,y,z) concerning colatitude
    mean(vecnorm([gradX(:,2),gradY(:,2),gradZ(:,2)], 2, 2));

    dS = vecnorm(cross([gradX(:,1),gradY(:,1),gradZ(:,1)],[gradX(:,2),gradY(:,2),gradZ(:,2)]), 2, 2);
    % scatter3(azi,col,dS) % Correct
    dS = dS ./ sin(col);
    % dS(sin(col) == 0) = 0;

    dO = 4*pi / size(azi,1); % dO is the average solid angle of the surface uniform mesh
    % dO = power(pi/size(col,1),2);
    S = sum(dS)*dO;
end