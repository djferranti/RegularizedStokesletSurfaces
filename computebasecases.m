function [t003, t001, se1m1, se2m1, sdm1, se1p1, se2p1, sdp1, ... 
    geometryData] = computebasecases(xField,Triangle,regularization)
%% COMPUTEBASECASES computes eight base cases for regularized stokeslet surfaces
% Parameters:
%   xField: 3 x M array of field points in 3D space
%   Triangle: struct containing data for a single triangle
%   regularization: regularization parameter for force spreading
%   geometryData: struct containing the geometry data needed for all the
%   other formulas
% Output:
%   t003: M x 1 vector representing T_{0,0,3} 
%   t001:  M x 1 vector representing T_{0,0,1} 
%   se1m1: M x 1 vector representing S^{e1}_{0,-1} 
%   se2m1: M x 1 vector representing S^{e2}_{0,-1} 
%   sdm1: M x 1 vector representing S^{d}_{0,-1} 
%   se1p1: M x 1 vector representing S^{e1}_{0,+1} 
%   se2p1: M x 1 vector representing S^{e2}_{0,+1}
%   sdp1: M x 1 vector representing S^{d}_{0,+1} 



M=size(xField,2); %M = number of field points
t003=zeros(M,1); t001 = t003; se1m1 = t003; se2m1 = t003; sdm1 = t003; 
se1p1 = t003; se2p1 = t003; sdp1 = t003;
bh = Triangle.bh; %for consistency with other T_{ijk}

%initialize these for geometryData
R0=zeros(M,1); R1=R0; R2=R0; 
x0DotV = zeros(M,1); x1DotW = x0DotV; x2DotD = x0DotV;
ell1 = 0; ell2 = 0; ell3 = 0;

for i = 1:3 %for every side of triangle

    %geometry data
    ya=Triangle.vertices(:,i); %vertex a of side i
    x0i=xField-ya; %vectors pointing from ya to field points
    nhatSidei=Triangle.normalstosides(:,i); %outward normal vector to side i
    sideiLength=Triangle.lengths(i); %length of side i
    vhati=Triangle.directions(:,i); %direction from yb to ya


    %repeat the vectors for dot product multiplication with x0
    nhatiRep=repmat(nhatSidei,1,M);
    vhatiRep=repmat(vhati,1,M);


    %more triangle geometry/field point data
    x0DotNi=(dot(x0i,nhatiRep))';
    x0DotVi=(dot(x0i,vhatiRep))';


    %gamma: squared distance to plane of triangle + squared regularization
    %only compute once - this is same for every iteration 
    %implementation note: numerical issues arose when computed on 
    % every iteration
    if i == 1
        r2Proj = dot(x0i,x0i)'-x0DotNi.^2-x0DotVi.^2; %squared dist to plane of triangle 
        r2Proj(r2Proj<eps)=0; %set to 0 if less than machine precision
        gamma = sqrt(r2Proj + regularization^2);

        %also compute x0DotW
        what = Triangle.directions(:,i+1); 
        whatRep = repmat(what,1,M); 
        x0DotW = (dot(x0i,whatRep))';

        %and save some other data relevant for i==1 case
        vhat = vhati;
        vDotW = dot(vhat, what);
        x0 = x0i;
    end


    t003 = t003 + computet003side(x0DotVi, x0DotNi, gamma, sideiLength);
    s0p1 = computes0p1(x0DotVi, x0DotNi, gamma, sideiLength);
    t001 = t001 + computet001side(s0p1, x0DotNi, sideiLength);
    s0m1 = computes0m1(x0DotVi, x0DotNi, gamma, sideiLength);

    se1p1 = se1p1 + (i==1) .*  s0p1; 
    se2p1 = se2p1 + (i==2) .*  s0p1; 
    sdp1  = sdp1  + (i==3) .*  s0p1;

    se1m1 = se1m1 + (i==1) .*  s0m1; 
    se2m1 = se2m1 + (i==2) .*  s0m1; 
    sdm1  = sdm1  + (i==3) .*  s0m1;

    R0 = R0 + (i==1) .* sqrt(dot(x0i,x0i)' + regularization.^2);
    R1 = R1 + (i==2) .* sqrt(dot(x0i,x0i)' + regularization.^2);
    R2 = R2 + (i==3) .* sqrt(dot(x0i,x0i)' + regularization.^2);

    x0DotV = x0DotV + (i==1) .* x0DotVi; 
    x1DotW = x1DotW + (i==2) .* x0DotVi; 
    x2DotD = x2DotD + (i==3) .* x0DotVi; 
    
    ell1 = ell1 + (i==1) .* sideiLength; 
    ell2 = ell2 + (i==2) .* sideiLength; 
    ell3 = ell3 + (i==3) .* sideiLength;


end %end for loop

%for t003 divide by gamma 
t003 = gamma.^(-1) .* t003; 
%t001 depends on t003 (see section on T001 in paper - there is a typo 
% epsilon^2 should be gamma^2 like it is here) 
t001 = t001 - gamma.^2 .* t003;

%for t003,t001 divide by BH since we multiply everything by BH in the end

t003 = bh.^(-1) .* t003; 
t001 = bh.^(-1) .* t001;

%put all the geometry data into a struct
geometryData = struct('x0', x0, 'vhat', vhat, 'what', what, 'x0DotV', ...
    x0DotV, 'x0DotW', x0DotW, 'x1DotW', x1DotW, 'x2DotD', x2DotD, ... 
    'R0', R0, 'R1', R1, 'R2', R2, 'vDotW', vDotW, 'ell1', ell1, ...
    'ell2', ell2, 'ell3', ell3, 'bh',bh);

end
