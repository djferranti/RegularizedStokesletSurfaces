function stokesletMatrix = ...
    assemblestokesletmatrix(xField,TriangleArray,numberTrianglePoints, ...
    regularization, mu) 
%% ASSEMBLESTOKESLETMATRIX assembles the regularized Stokeslet surface matrix. 
% Parameters:
%   xField: 3 x M array of field points 
%   TriangleArray: 1 x Q struct array where Q is the number of triangular
%   faces. 
%   numberTrianglePoints: number of unique points that make up
%   triangulation
%   regularization: blob parameter
%   mu: viscosity parameter
% Output:
%   A: 3M x 3V regularized Stokeslet surface matrix which is related to the
%   velocity U at the field points xf by U = A*F where F are the forces at 
%   the V distinct vertices of the Q triangular faces.

numberFaces = size(TriangleArray,2); 
numberFieldPoints = size(xField,2);

stokesletMatrix = zeros(3*numberFieldPoints, 3*numberTrianglePoints);

for q = 1: numberFaces
    %triangle q
    Triangle = TriangleArray(q);
    bh = Triangle.bh;


    %compute the base cases
    [t003,t001,se1m1,se2m1,sdm1,se1p1,se2p1,sdp1, geometryData] = ...
        computebasecases(xField, Triangle, regularization);

    %compute the other T_{mnq} recursively
    [t101, t011] = tqequals1(se1m1, se2m1, sdm1, t001, geometryData);
    [t103, t013, t203, t023, t113, t303, t033, t213, t123] = ...
        tqequals3(se1p1, se2p1, sdp1, se1m1, se2m1, sdm1, t001, t003, ...
        t101, t011, geometryData);

    %compute the block columns for the stokeslet matrix
    [b0,b1,b2] = computeblocks(geometryData, t001, t101, t011, ...
        t003, t103, t013, t203, t023, t113, t303, t033, t213, t123, regularization);

    %put the blocks into relevant block columns corresponding to indices
    indices = Triangle.indices;
    index0 = indices(1); index1 = indices(2); index2 = indices(3);

    addBlockColumn(b0,index0);
    addBlockColumn(b1,index1);
    addBlockColumn(b2,index2);
    
end

    function addBlockColumn(block, index)
        stokesletMatrix(:, 3 * (index - 1) + 1: 3 * (index - 1) + 1 + 2) = ...
            stokesletMatrix(:, 3 * (index - 1) + 1 : 3 * (index - 1) + 1 + 2) + ...
            block .* (bh ./ (8 * pi * mu));
    end

end