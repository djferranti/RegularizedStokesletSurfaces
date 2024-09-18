function [a101, a011, a201, a021, b101, b011, b201, b021] = ...
    abmnq(a001, se10p1, se20p1, sd0p1, ... 
    se10m1, se20m1, sd0m1, geometryData)
%% ABMNQ evaluates A101, A011, A201, A021, and their B counterparts
% Parameters:
%   a001: M x 1 vector storing A_{0,0,1}
%   se10p1, se20p1, sd0p1: M x 1 vectors storing the integral evaluation for
%   S_{0,1} along the respective sides in directions e_1, e_2, d.
%   se10m1, se20m1, sd0m1: M x 1 vectors storing the integral evaluation for
%   S_{0,-1} along the respective sides in directions e_1, e_2, d.
%   geometryData: struct containing data for x0DotV, x1DotW, x2DotD, L1,
%   L2, L3
%   Output:
%   a001, a011, a201, a021, b101, b011, b201, b021: M x 1 vectors storing
%   A_{1,0,1], A_{0,1,1}, A_{2,0,1}, A_{0,2,1}, B_{1,0,1}, B_{0,1,1], 
%   B_{2,0,1}, B_{0,2,1}


%unpackage geometryData (M x 1 vectors)
x0DotV = geometryData.x0DotV;
% x0DotW = geometryData.x0DotW;
% vDotW = geometryData.vDotW;
x1DotW = geometryData.x1DotW; 
x2DotD = geometryData.x2DotD;
R0 = geometryData.R0; 
R1 = geometryData.R1; 
R2 = geometryData.R2;
%the rest are scalars
ell1 = geometryData.ell1; 
ell2 = geometryData.ell2;
ell3 = geometryData.ell3;

%compute the line integral recursions from formula (2.24) 
q = 1; %corresponds to q from formula, and easier on eyes

se1p1p1 = 1 ./ (ell1 .^ 2) .* ( -1 ./ (q-2) .* (R1 - R0) - ...
    ell1 .* x0DotV .* se10p1 );
se1p2p1 = 1 ./ (ell1 .^ 2) .* ( -1 ./ (q-2) .* R1 + 1 ./ (q-2) .* ...
    se10m1 - ell1 .* x0DotV .* se1p1p1 );

se2p1p1 = 1 ./ (ell2 .^ 2) .* ( -1 ./ (q-2) .* (R2 - R1) - ...
    ell2 .* x1DotW .* se20p1 );
se2p2p1 = 1 ./ (ell2 .^ 2) .* ( -1 ./ (q-2) .* R2 + 1 ./ (q-2) .* ...
    se20m1 - ell2 .* x1DotW .* se2p1p1 );

sdp1p1 = 1 ./ (ell3 .^ 2) .* ( -1 ./ (q-2) .* (R0 - R2) - ...
    ell3 .* x2DotD .* sd0p1 );
sdp2p1 = 1 ./ (ell3 .^ 2) .* ( -1 ./ (q-2) .* R0 + 1 ./ (q-2) .* ...
    sd0m1 - ell3 .* x2DotD .* sdp1p1 );

%compute the Amnq, Bmnq from formulas (2.22), (2.23) 
%formula (2.22) in paper should have had a (m+n) choose k coefficient
a101 = a001 + sdp1p1; 
a201 = a001 + 2*sdp1p1 - sdp2p1; 
a011 = se2p1p1 - sd0p1 + sdp1p1; 
a021 = se2p2p1 - sd0p1 + 2*sdp1p1 - sdp2p1; 

%this case appears to be unnecessary
% a111 = se2p1p1 - sd0p1 + sdp1p1 - sdp2p1;

%n=0 formula part of (2.23)
b101 = -se1p1p1 + sd0p1 - sdp1p1;
b201 = -se1p2p1 + sd0p1  - 2*sdp1p1 + sdp2p1;

%n>0 formula part of (2.23) 
b011 = sd0p1 - sdp1p1;
b021 = sd0p1 - 2*sdp1p1 + sdp2p1; 

%this case appears also to be unnecessary
%b111 = b021; 



% debug stuff
% BH = geometryData.BH;
% nsdp2p1 = nIntegrateLineAlpha2Rm1(x2DotD,ell3,R2);
% nsdp2p1 = nsdp2p1./ell3; 
% [sdp2p1, nsdp2p1]
% 
% nse1p2p1 = nIntegrateLineAlpha2Rm1(x0DotV,ell1,R0);
% nse1p2p1 = nse1p2p1./ell1; 
% [se1p2p1, nse1p2p1];

% na001 = nIntegrateA001(x0DotV,x0DotW,vDotW,ell1,ell2,R0); 
% na001 = na001./BH;
% [a001, na001]

% m=1;
% na101 = nIntegrateAm01(m,x0DotV,x0DotW,vDotW,ell1,ell2,R0);
% na101 = na101./BH;
% [a101, na101]
% 
% m=2;
% na201 = nIntegrateAm01(m,x0DotV,x0DotW,vDotW,ell1,ell2,R0);
% na201 = na201./BH;
% [a201, na201]

function out=nIntegrateA001(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    Rsq=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    am01= @(a,b) - ( (x0DotV(i) + ell1.*a).*ell1 ...
        + ell1.*ell2.*b.*vDotW ) .* Rsq(a,b).^(-3/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(am01,0,1,0,bmax);
end%compare with numerical integration
end

function out=nIntegrateAm01(m,x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    Rsq=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    am01= @(a,b) m .* a.^(m-1) .* Rsq(a,b).^(-1/2) - a.^m .* ( (x0DotV(i) + ell1.*a).*ell1 ...
        + ell1.*ell2.*b.*vDotW ) .* Rsq(a,b).^(-3/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(am01,0,1,0,bmax);
end%compare with numerical integration
end

function out=nIntegrateLineAlphaRm1(x0DotV,ell1,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    Rsq=@(a)(x0DotV(i)+a.*ell1).^2 + R0(i).^2 - (x0DotV(i)).^2;
    Rm1=@(a) a.* Rsq(a).^(-1/2);
    out(i) = ell1.*integral(Rm1,0,1);
end
end


function out=nIntegrateLineAlpha2Rm1(x0DotV,ell1,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    Rsq=@(a)(x0DotV(i)+a.*ell1).^2 + R0(i).^2 - (x0DotV(i)).^2;
    Rm1=@(a) a.^2 .* Rsq(a).^(-1/2);
    out(i) = ell1.*integral(Rm1,0,1);
end
end

end