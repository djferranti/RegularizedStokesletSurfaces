function s0p1Side = computes0p1(x0DotV, x0DotN, gamma, sideLength)
%% COMPUTES0p1 evaluates line integral S_{0,1} for side of a triangle
% Parameters:
%   x0DotV: M x 1 vector representing x(0) dot v for every field point, where 
%   v is direction from endpoint to start point of side
%   x0DotN: M x 1 vector representing x0 dot n for every field point, where
%   n is unit outward normal to 
%   gamma: M x 1 vector representing regularized distance to plane of
%   triangle for every field point
%   sideLength: length of side (scalar)
%
% Output:
%   s0p1Side: M x 1 vector representing the line integral S_{0,1}

    
s0p1Side = formula226;


%% helper functions

    function out = formula226 %formula (2.26) in paper
        x1DotV = x0DotV+sideLength; %x0DotV+L
        R1 = sqrt(x1DotV.^2 + x0DotN.^2 + gamma.^2);
        R0 = sqrt(x0DotV.^2 + x0DotN.^2 + gamma.^2);

        out = 1 / sideLength .* ( atanh(x1DotV./R1) - atanh(x0DotV./R0) );
    end
end
