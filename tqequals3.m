function [t103, t013, t203, t023, t113, t303, t033, t213, t123] = ... 
    tqequals3(se10p1, se20p1, sd0p1, se10m1, se20m1, sd0m1, ... 
    t001,t003, t101, t011, geometryData) 
%% TQEQUALS3 evaluates T103, T013, T203, T023, T113, T303, T033, T213, T123
% Parameters:
%   se10p1, se20p1, sd0p1: M x 1 vectors storing the integral evaluation for
%   S_{0,1} along the respective sides in directions e_1, e_2, d.
%   se10m1, se20m1, sd0m1: M x 1 vectors storing the integral evaluation for
%   S_{0,-1} along the respective sides in directions e_1, e_2, d.
%   t001, t003, t101, t011: M x 1 vector storing the integral evaluation for 
%   T_{0,0,1}, T_{0,0,3}, T_{1,0,1}, T_{0,1,1} 
%   geometryData: struct containing data for x0DotV, x0DotW, vDotW, L1, L2
%
% Output:
%   t103, t013, t203, t023, t113, t303, t033, t213, t123: M x 1 vectors
%   storing the integral evaluations of T_{1,0,3}, T_{0,1,3}, T_{2,0,3},
%   T_{0,2,3}, T_{1,1,3}, T_{3,0,3}, T_{0,3,3}, T_{2,1,3}, T_{1,2,3}

%initialize data structures to be output
M = size(t001,1); 

%formulas (2.22) and (2.23) 
a001 = se20p1 - sd0p1;
b001 = -se10p1 + sd0p1; 

[a101, a011, a201, a021, b101, b011, b201, b021] = ... 
    abmnq(a001, se10p1, se20p1, sd0p1, ...
    se10m1, se20m1, sd0m1, geometryData);

%these are useful when m=0 or n=0, allows for only one recursion formula
tmm1qm2Dummy = zeros(M,1); tnm1qm2Dummy = zeros(M,1);

%use recursion to move to stage m+n=1
m = 0; n=0; q=3; 
t103 = tmp1recursion(m,n,q,a001,b001,tmm1qm2Dummy,tnm1qm2Dummy,t003, ... 
    geometryData);
t013 = tnp1recursion(m,n,q,a001,b001,tmm1qm2Dummy,tnm1qm2Dummy,t003, ... 
    geometryData);


%use recursion to move to stage m+n=2
m=1; n=0; q=3; 
t203 = tmp1recursion(m,n,q,a101,b101,t001,tnm1qm2Dummy,t103, ... 
    geometryData); 

t113 = tnp1recursion(m,n,q,a101,b101,t001,tnm1qm2Dummy,t103, ...
    geometryData);

m=0; n=1; q=3; 
t023 = tnp1recursion(m,n,q,a011,b011,tmm1qm2Dummy,t001,t013, ... 
    geometryData);

%use recursion to move to stage m+n=3
m=2; n=0; q=3; 
t303 = tmp1recursion(m,n,q,a201,b201,t101,tnm1qm2Dummy,t203, ... 
    geometryData); 

t213 = tnp1recursion(m,n,q,a201,b201,t101,tnm1qm2Dummy,t203, ... 
    geometryData); 

m=0; n=2; q=3; 
t033 = tnp1recursion(m,n,q,a021,b021,tmm1qm2Dummy,t011,t023, ... 
    geometryData); 

t123 = tmp1recursion(m,n,q,a021,b021,tmm1qm2Dummy,t011,t023, ... 
    geometryData); 

end