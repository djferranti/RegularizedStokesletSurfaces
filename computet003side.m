function t003Side = computet003side(x0DotV, x0DotN, gamma, sideLength)
%% COMPUTET003SIDE evaluates contribution from side in contour integral equivalent to T003
% Parameters:
%   x0DotV: M x 1 vector representing x0 dot v for every field point, where 
%   v is direction from endpoint to start point of side
%   x0DotN: M x 1 vector representing x0 dot n for every field point, where
%   n is unit outward normal to 
%   gamma: M x 1 vector representing regularized distance to plane of
%   triangle for every field point
%   sideLength: length of side (scalar)
%
% Output:
%   t003Side: M x 1 vector representing contribution from side in the 
%   contour integral equivalent to T003

%from analytic integration formula - for details, see paper
p=x0DotV./sideLength; 
q=sqrt((x0DotN./sideLength).^2 ...
    +(gamma.^2)./(sideLength).^2); 
x0DotNneq0 = abs(x0DotN)>eps ; 
isInPInterval = (-1 < p) & (p < 0); 
condSet1 = x0DotNneq0 & isInPInterval;
condSet2 = x0DotNneq0 & ~isInPInterval;


%this is different from paper - the paper should have contained a
%multiplicative factor of L/H which turns it into this.
integralCoefficient=(-1) .* x0DotN ./ sideLength;


%the following are used several times in the helper functions below
[r1,r2] = r1r2; 
ellq = sideLength.*q;
onePlus = 1 + gamma./ellq;
oneMinus = 1 - gamma./ellq;
sqrtOnePlusMinus = sqrt(onePlus ./ oneMinus);
sqrtOneMinusPlus = sqrt(oneMinus ./ onePlus);

%1st/2nd formula only modifies t003 when condSet1/2 true, x0DotN neq 0
t003Side = condSet1.*integralCoefficient.*formula230 + ...
    condSet2.*integralCoefficient.*formula229;
t003Side(abs(x0DotN) < eps) = 0; 

%another special case is when 1 - gamma./ellq = 0 
%in the limit as this goes to zero, the resulting expression is 0
t003Side(abs(oneMinus) < eps) = 0;

%% helper functions
function [r1,r2]=r1r2
        phi1 = asec( sqrt(p.^2./q.^2 + 1) )./2;
        phi2 = asec( sqrt((1+p).^2./q.^2 + 1) )./2;
        r1 = tan(phi1);
        r2 = tan(phi2);
end

function intCase1 = formula230 %formula (2.30) in paper

%integral for p between -1 and 0
intCase1 = 2 ./ (q.*(onePlus)) .* sqrtOnePlusMinus .* ( ...
    atan(r1 .* sqrtOneMinusPlus) + atan(r2 .* sqrtOneMinusPlus));
end

function intCase2 = formula229 %formula (2.29) in paper

%integral for p <= -1 or >= 0
sgnfac = sign(1+p);
sgnfac(abs(sgnfac) < eps) = -1;
intCase2 = sgnfac .* 2 ./ (q.*(onePlus)) .* sqrtOnePlusMinus .* ...
    (atan(r2.* sqrtOneMinusPlus) - atan(r1 .* sqrtOneMinusPlus));
end
end