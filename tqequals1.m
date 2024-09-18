function [t101, t011] = tqequals1(se1m1, se2m1, sdm1, t001, geometryData) 
%% TQEQUALS1 outputs T101, T011 from recursion formula
% Parameters:
%   se1m1, se2m1, sdm1: M x 1 vectors storing the integral evaluation
%   of S_{0,1} along the respective sides in directions e_1, e_2, d.
%   t001: M x 1 vector storing the integral evaluation of T_{0,0,1} 
%   geometryData: struct containing data for x0DotV, x0DotW, vDotW, L1, L2
%
% Output:
%   t101: M x 1 vector storing the integral evaluation of T_{1,0,1}
%   t011: M x 1 vector storing the integral evaluation of T_{0,1,1}

%initialize data structures to be output
M = size(t001,1); 

%formulas (2.22) and (2.23) 
a00m1 = se2m1 - sdm1;
b00m1 = -se1m1 + sdm1; 

%these are irrelevant in this particular recursion
tm10m2 = zeros(M,1); t0m1m2 = zeros(M,1); 

%recursion indices
m=0; n=0; q=1;

t101 = tmp1recursion(m,n,q,a00m1,b00m1,tm10m2,t0m1m2,t001,geometryData);
t011 = tnp1recursion(m,n,q,a00m1,b00m1,tm10m2,t0m1m2,t001,geometryData);

end