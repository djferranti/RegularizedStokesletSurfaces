function [TriangleArray, xyzPts, faces] =  triangulatesphereicos(factor, radius, varargin)
%% TRIANGULESPHEREICOS creates array of triangle structs from icosahedral triangulation of sphere
% Parameters:
%   factor: number of subdivisions of the icosahedron used.
%   radius: radius of sphere
%   varargin: index start
%
% Output:
%   TriangleArray: Struct array which contains triangle structs with
%   relevant data for regularized Stokeslet surfaces
%   xyzPts: 3 x V array of points where V is the number of distinct
%   vertices in the triangulation of the sphere.
%   faces: 3 x N array of indices where N is the number of triangular
%   faces.
%
% Attribution: Uses John Burkardt's sphere_delaunay code distributed
% under the GNU LGPL license. Details on the implementation can be found at
% the website
% https://people.math.sc.edu/Burkardt/m_src/sphere_delaunay/sphere_delaunay.html.

if nargin == 3
    indexStart = varargin{1};
else 
    indexStart = 0;
end

%add the relative path to Burkhardt code
% addpath(".\sphere_delaunay\")

numberPts = sphere_grid_icos_size(factor);
xyzPts = radius .* sphere_gridpoints_icos2(factor,numberPts);


%compute the Delaunay triangulation of the points.
[~, faces] = sphere_delaunay(numberPts, xyzPts);
faces = faces + indexStart;
TriangleArray = createtrianglestruct(xyzPts,faces);

    function TriangleArray=createtrianglestruct(xyzPts,faces)
        nfaces=size(faces,2);
        for ff=1:nfaces
            indices = faces(:,ff)';
            pt1 = xyzPts(:,indices(1) - indexStart); 
            pt2 = xyzPts(:,indices(2) - indexStart);
            pt3 = xyzPts(:,indices(3) - indexStart);
            normal=cross(pt1 - pt2, pt2 - pt3);

            %check if normal points inside to center of sphere. 
            % if so, change pt1,pt3.
            if dot(normal, pt1) < 0
                pt1Prev = pt1; 
                pt1 = pt3;
                pt3 = pt1Prev;
                indices = [indices(3), indices(2), indices(1)];
                faces(:,ff) = indices';
            end

            vertices=[pt1,pt2,pt3];
            directions=[pt1-pt2,pt2-pt3,pt3-pt1];
            lengths=vecnorm(directions);
            directions=directions./lengths;
            normaltoplane=cross(directions(:,1),directions(:,2));
            normaltoplane=normaltoplane./vecnorm(normaltoplane);
            normalstosides=[cross(normaltoplane,directions(:,1)), ...
                cross(normaltoplane,directions(:,2)),...
                cross(normaltoplane,directions(:,3))];
            ht1=lengths(2).*dot(directions(:,2),normalstosides(:,1));
            bh=lengths(1).*ht1;
            heights=[ht1,bh./lengths(2),bh./lengths(3)];

            Triangle=struct('vertices',vertices, 'indices', indices, ...
                'directions',directions, 'normaltoplane',normaltoplane, ...
                'normalstosides',normalstosides,'lengths',lengths, ...
                'heights',heights,'bh',bh);
            TriangleArray(ff) = Triangle; 
        end
    end
end
