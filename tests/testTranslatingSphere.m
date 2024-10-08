%testTranslatingSphere 
%the "forward problem": given force density over sphere, evaluate velocity
addpath(genpath('./..'))
regularization = 1e-6; 
mu = 1;
%use "factor" number of divisions of the icosahedron to triangulate sphere
factors = [3, 6, 8, 12];
% factors = [3, 6, 8, 12, 24]; %this one may take a while on a desktop
radius = 1; 
ell2Errors = zeros(size(factors)); 
relativeErrors = zeros(size(factors));
averageDistance = zeros(size(factors));

for i = 1 : length(factors)
    %triangulate sphere and put triangle struct in
    % TriangleArray - a 1 x numberOfTriangles struct array.
    %see triangleulatesphereicos function to understand the data needed 
    % in a Triangle struct
    ithFactor = factors(i);
    [TriangleArray, points, faces] = triangulatesphereicos(ithFactor, ...
        radius);
    
    % if you want to see sphere
    % plotSphere(TriangleArray, points, faces)

    xField = points; %use triangle points as field points
    numberTrianglePoints = size(points,2);
    disp(['number of triangles = ' num2str(size(TriangleArray,2))])
    disp(['number of DOF = ' num2str(3 * numberTrianglePoints)])
    
    %assemble Stokeslet matrix
    A = assemblestokesletmatrix(xField,TriangleArray,numberTrianglePoints, ...
        regularization, mu); 

    disp(['condition number of A = ' num2str(cond(A))])


    %impose constant traction on sphere in x-direction
    unitX = [1,0,0]';
    F = 3 * mu / (2 * radius) * unitX;
    F = repmat(F, 1, numberTrianglePoints);
    %use (:) to write F as one long vector
    F = F(:);

    %multiply stokeslet matrix by F to get the velocity at each point
    uField = A*F; 

    %reshape into 3 x numberTrianglePoints array
    uField = reshape(uField,3,numberTrianglePoints);

    %the traction imposes a unit velocity in the theoretical case
    %check that the discretized version is approximating this

    uTheory = repmat(unitX, 1, numberTrianglePoints); 
    error = uField - uTheory; 
    ell2Errors(i) = sqrt( sum( sum( error .* error ) ) / numberTrianglePoints ) 

end

function plotSphere(TriangleArray, points, faces)
    hold off
    patch('Faces',faces','Vertices',points', ...
    'FaceColor','c', 'FaceAlpha', 0);
    hold on
    for jj = 1:size(TriangleArray,2) 
        Triangle = TriangleArray(jj);
        centroid = sum(Triangle.vertices, 2) ./ 3; 
        xc = centroid(1); yc = centroid(2); zc=centroid(3); 
        normalToPlane = Triangle.normaltoplane;
        nx = normalToPlane(1); ny =normalToPlane(2); nz = normalToPlane(3);
        quiver3(xc,yc,zc,nx,ny,nz,0.5,'k');
    end
    axis equal
    hold off
end