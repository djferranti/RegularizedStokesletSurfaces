function velos = ...
    evaluatevelocity(xField,xyzTriangle,forcesTriangle,TriangleArray, ...
    regularization, mu) 
%% EVALUATEVELOCITY evaluates the velocity due to the force density over triangulated surface.
% Parameters:
%   xField: 3 x M array of field points 
%   xyzTriangle: 3 x N array of triangle points
%   forcesTriangle: 3 x N array of forces at corresponding triangle points
%   TriangleArray: 1 x Q struct array where Q is the number of triangular
%   faces. 
%   triangulation
%   regularization: blob parameter
%   mu: viscosity parameter
% Output:
%   velos: 3 x M array of velocities at the M field points 

numberFaces = size(TriangleArray,2); 
numberFieldPoints = size(xField,2);

velos = zeros(3, numberFieldPoints);

for q = 1: numberFaces
    %triangle q
    Triangle = TriangleArray(q);
    bh = Triangle.bh; 
    f0 = forcesTriangle(:, Triangle.indices(1)); 
    f1 = forcesTriangle(:, Triangle.indices(2)); 
    f2 = forcesTriangle(:, Triangle.indices(3)); 
    fa = f1 - f0; fb = f2 - f1;

    y0 = xyzTriangle(:, Triangle.indices(1));

    [p00,p10,p01,p20,p11,p02,p30,p21,p12,p03] = computepcoeffs(xField, ...
        f0, f1, f2, y0, Triangle, regularization);

    %compute the base cases
    [t003,t001,se1m1,se2m1,sdm1,se1p1,se2p1,sdp1, geometryData] = ...
        computebasecases(xField, xyzTriangle, Triangle, regularization);

    %compute the other T_{mnq} recursively
    [t101, t011] = tqequals1(se1m1, se2m1, sdm1, t001, geometryData);
    [t103, t013, t203, t023, t113, t303, t033, t213, t123] = ...
        tqequals3(se1p1, se2p1, sdp1, se1m1, se2m1, sdm1, t001, t003, ...
        t101, t011, geometryData); 

    f0Rep = repmat(f0,1,numberFieldPoints);
    faRep = repmat(fa,1,numberFieldPoints); 
    fbRep = repmat(fb,1,numberFieldPoints);

    velos = velos + bh ./ (8 * pi * mu) .* (...
        f0Rep .* t001' + ... 
        p00 .* t003' + faRep .* t101' + ...
        p10 .* t103' + fbRep .* t011' + ... 
        p01 .* t013' + p20 .* t203' + ...
        p11 .* t113' + p02 .* t023' + ...
        p30 .* t303' + p21 .* t213' + ...
        p12 .* t123' + p03 .* t033' ...
        ) ;

end
