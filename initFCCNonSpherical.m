solidVolumeFraction = 0.2;
xFeret = 0.1;
yFeret = 0.2;
zFeret = 0.5;
NUC = 3;
volume = pi*xFeret*yFeret*zFeret/6;

% Defining FCC unit cell 
fccbase = [...    
       0    0    0  ;...
       0    0    1  ;...
       0    1    0  ;...
       0    1    1  ;...
       1    0    0  ;...
       1    0    1  ;...
       1    1    0  ;...
       1    1    1  ;...
       0.5  0.5  0  ;...
       0.5  0.5  1  ;...
       0.5  0    0.5;...
       0.5  1    0.5;...
       0    0.5  0.5;...
       1    0.5  0.5];

[latx,laty,latz] = ndgrid(0:NUC-1,0:NUC-1,0:NUC-1);
latxyz = repmat([latx(:),laty(:),latz(:)],14,1);
latxyz = latxyz + reshape(permute(repmat(fccbase,[1,1,numel(latx)]),[3 1 2]),[],3);
latxyz = unique(latxyz,'rows');

Lx = power(size(latxyz,1)*volume*power(xFeret,2)/(solidVolumeFraction*yFeret*zFeret), 1/3);
alpha = (Lx*0.99 - xFeret) / NUC;
latxyz(:,1) = alpha*latxyz(:,1) + xFeret/2;

Ly = power(size(latxyz,1)*volume*power(yFeret,2)/(solidVolumeFraction*xFeret*zFeret), 1/3);
alpha = (Ly*0.99 - yFeret) / NUC;
latxyz(:,2) = alpha*latxyz(:,2) + yFeret/2;

Lz = power(size(latxyz,1)*volume*power(zFeret,2)/(solidVolumeFraction*yFeret*xFeret), 1/3);
alpha = (Lz*0.99 - zFeret) / NUC;
latxyz(:,3) = alpha*latxyz(:,3) + zFeret/2;

plotSpheres(latxyz(:,1),latxyz(:,2),latxyz(:,3),xFeret,yFeret,zFeret)

function S = plotSpheres(x,y,z,xFeret,yFeret,zFeret)
    s = size(x,1);
    [unitX, unitY, unitZ] = ellipsoid(0.0, 0.0, 0.0, xFeret/2, yFeret/2, zFeret/2);
    for i = 1 : s
      X = unitX + x(i);
      Y = unitY + y(i);
      Z = unitZ + z(i);
      surf(X,Y,Z);
      hold on
    end  
end