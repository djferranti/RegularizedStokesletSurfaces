function [p00,p10,p01,p20,p11,p02,p30,p21,p12,p03] = computepcoeffs(xField, ...
        f0, f1, f2, y0, Triangle, regularization) 
%% COMPUTEPCOEFFS computes the vector coefficients for velocity evaluation due to force over triangle patch.
% Parameters:
%   xField: 3 x M array of field points 
%   f0,f1,f2: 3 x 1 vectors corresponding to forces at triangle points y0,
%   y1, y2
%   y0: the "starting" vertex of the triangle
%   Triangle: struct with geometric data for triangle
%   regularization: blob parameter
% Output:
%   p00,p10,p01,p20,p11,p02,p30,p21,p12,p03: 3 x M array of coefficients 
%   required for velocity evaluation at field points

M = size(xField,2);

fa = f1 - f0; fb = f2 - f1; 

vhat = Triangle.directions(:,1); 
what = Triangle.directions(:,2); 
ell1 = Triangle.lengths(1);
ell2 = Triangle.lengths(2);

vhatRep = repmat(vhat,1,M);
whatRep = repmat(what,1,M);

f0DotV = dot(f0,vhat);
f0DotW = dot(f0,what);
faDotV = dot(fa,vhat); 
faDotW = dot(fa,what); 
fbDotV = dot(fb,vhat);
fbDotW = dot(fb,what); 

x0 = xField - y0;

f0DotX0 = repmat ( sum (repmat(f0,1,M) .* x0, 1), 3, 1);
faDotX0 = repmat ( sum (repmat(fa,1,M) .* x0, 1), 3, 1);
fbDotX0 = repmat ( sum (repmat(fb,1,M) .* x0, 1), 3, 1);

reg2 = regularization .^ 2;

p00 = reg2 .* f0 + (f0DotX0) .* x0; 
p10 = reg2 .* fa + (ell1.*f0DotV + faDotX0) .* x0 + (ell1 .* f0DotX0) .* vhatRep;
p01 = reg2 .* fb + (ell2.*f0DotW + fbDotX0) .* x0 + (ell2 .* f0DotX0) .* whatRep; 

p20 = (ell1 .* faDotV) .* x0 + (ell1 .^2 .* f0DotV + ell1 .* faDotX0 ) .* vhatRep; 
p02 = (ell2 .* fbDotW) .* x0 + (ell2 .^2 .* f0DotW + ell2 .* fbDotX0 ) .* whatRep;
p11 = (ell1 .* fbDotV + ell2 .* faDotW) .* x0 + (ell1 .* ell2 .* f0DotW  + ell1 .* fbDotX0) .* vhatRep ...
    + (ell1 .* ell2 .* f0DotV + ell2 .* faDotX0) .* whatRep;

p30 = (ell1 .^2 .* faDotV) .* vhatRep; 
p03 = (ell2 .^2 .* fbDotW) .* whatRep;
p21 = (ell1 .* ell2 .* faDotW + ell1 .^2 .* fbDotV) .* vhatRep + (ell1 .* ell2 .* faDotV) .* whatRep;
p12 = (ell1 .* ell2 .* fbDotV + ell2 .^2 .* faDotW) .* whatRep + (ell1 .* ell2 .* fbDotW) .* vhatRep;


end