%testDragSphere
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


    %impose unit velocity of sphere in x direction
    unitX = [1,0,0]';
    uField = repmat(unitX, 1, numberTrianglePoints);


    %use (:) to write uField as [1,0,0,1,0,0,...]'
    uField = uField(:);

    %solve for forces F then reshape into 3 x numberTrianglePoints array
    F = A \ uField;
    F = reshape(F,3,numberTrianglePoints);

    %organize forces F0,F1,F2 into 3 x numberOfTriangles array, so that column j
    %of F0 array corresponds to the force at y0 of the jth triangle, etc.
    F0=F(:,faces(1,:)); F1=F(:,faces(2,:)); F2=F(:,faces(3,:));

    %output base times height of each triangle to calculate drag
    bh = [TriangleArray.bh]; %this is how you output comma separated lists!
    totalDrag=sum((F0+F1+F2)./3.*bh./2,2); %drag formula from paper

    error = totalDrag - [6 * pi * mu * radius, 0, 0]'; 
    ell2Errors(i) = sqrt( sum( error .* error ) )
    relativeErrors(i) = abs(1 - totalDrag(1) ./ (6 * pi* mu * radius) ) ...
        .* 100
    averageDistance(i) = sqrt(mean(bh));

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