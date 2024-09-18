function s0m1Side = computes0m1(x0DotV, x0DotN, gamma, sideLength)
%% COMPUTES0m1 computes line integral S_{0,-1} for side of a triangle
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
%   s0m1Side: M x 1 vector representing the line integral S_{0,-1}

    
s0m1Side = formula225;


%% helper functions

    function out = formula225 %formula (2.25) in paper
        x1DotV = x0DotV+sideLength; %x0DotV+L
        R1 = sqrt(x1DotV.^2 + x0DotN.^2 + gamma.^2);
        R0 = sqrt(x0DotV.^2 + x0DotN.^2 + gamma.^2);

        out = 1 / (2*sideLength) .* ( ( x1DotV.*R1+(x0DotN.^2 + gamma.^2)... 
            .*log(x1DotV + R1) ) -  ( x0DotV.*R0 + (x0DotN.^2 + gamma.^2).* ...
            log(x0DotV+R0) ) );

    end
end
