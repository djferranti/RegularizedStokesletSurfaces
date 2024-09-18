function t001Side = computet001side(s0p1, x0DotN, sideLength)
%% COMPUTET001SIDE computes contribution from side in contour integral part of T001
% Parameters:
%   s0p1: M x 1 vector representing the line integral S_{0,1} for each
%   field point
%   x0DotN: M x 1 vector representing x0 dot n for every field point, where
%   n is unit outward normal to side 
%   sideLength: length of side (scalar)
%
% Output:
%   t001Side: M x 1 vector representing contribution from side in the 
%   contour integral part of T001
    
%multiply by sideLength since the formula for s0p1 has a 1/sideLength
t001Side = sideLength .* ( (-1) .* (x0DotN .* s0p1 ) );
%if |x0DotN|<eps, then zero out the contribution from this side
t001Side(abs(x0DotN) < eps) = 0;