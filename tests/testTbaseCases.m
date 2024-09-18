%test Tmnq evaluations
addpath(genpath('./..'))

%create 100 evaluation points in unit cube
N=100;
rng(194)
xf=rand(3,N);
%add vertices of triangle as well  
% (these cases demonstrate the difficulty of the near singular integrals
% using default numerical integration techniques)
xf=[xf,Triangle.vertices];
N=N+3;
reg=10^(-6);

%exact formulas
[t003,t001,se1m1,se2m1,sdm1,se1p1,se2p1,sdp1, geometryData] = computebasecases(xf, Triangle, reg);
[t101, t011] = tqequals1(se1m1, se2m1, sdm1, t001, geometryData); 
[t103, t013, t203, t023, t113, t303, t033, t213, t123] = ...
    tqequals3(se1p1, se2p1, sdp1, se1m1, se2m1, sdm1, t001, t003, ...
    t101, t011, geometryData);

%numerical integration
%data needed for comparison (a lot of this is also in geometryData)
y0=Triangle.vertices(:,1); y1 = Triangle.vertices(:,2); 
y2 = Triangle.vertices(:,3);
x0=xf-y0; %vectors pointing from ya to field points
x1=xf-y1; x2=xf-y2;
R0=(sqrt(dot(x0,x0)+reg.^2))'; %regularized distance
R1=(sqrt(dot(x1,x1)+reg.^2))'; 
R2=(sqrt(dot(x2,x2)+reg.^2))'; %regularized distance
nhatSide=Triangle.normalstosides(:,1); %outward normal vector to side i
ell1=Triangle.lengths(1); %length of side 1
ell2=Triangle.lengths(2); %length of side 2 
ell3=Triangle.lengths(3); %length of side 2 
vhat=Triangle.directions(:,1); %direction of side 1 (y1 to y0)
what=Triangle.directions(:,2); %direction of side 2 (y2 to y1)
dhat=Triangle.directions(:,3); %direction of side 3 (y0 to y2)
vDotW=dot(vhat,what); 
bh = Triangle.bh;

%repeat the vectors for dot product multiplication with x0
vhatRep=repmat(vhat,1,N);
whatRep=repmat(what,1,N);
dhatRep=repmat(dhat,1,N);

%more triangle geometry/field point data + transpose arrays
x0DotW=(dot(x0,whatRep))';
x0DotV=(dot(x0,vhatRep))';
x1DotW=(dot(x1,whatRep))';
x2DotD=(dot(x2,dhatRep))';

%gamma^2 = squared distance to plane of triangle + squared regularization
gamma = sqrt((dot(x0,x0))'-x0DotW.^2-x0DotV.^2 +reg^2);

%numerical integration 
nt003=nIntegrateRm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0);
nt001=nIntegrateRm1(x0DotV,x0DotW,vDotW,ell1,ell2,R0);
nt101=nIntegrateAlphaRm1(x0DotV,x0DotW,vDotW,ell1,ell2,R0); 
nt011=nIntegrateBetaRm1(x0DotV,x0DotW,vDotW,ell1,ell2,R0); 
nt103=nIntegrateAlphaRm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0); 
nt013=nIntegrateBetaRm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0);
nt203=nIntegrateAlpha2Rm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0); 
nt023=nIntegrateBeta2Rm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0); 
nt113=nIntegrateAlphaBetaRm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0); 
nt303=nIntegrateAlpha3Rm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0); 
nt033=nIntegrateBeta3Rm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0); 
nt213=nIntegrateAlpha2BetaRm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0); 
nt123=nIntegrateAlphaBeta2Rm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0); 

%for comparison need to divide by BH
nt003 = nt003./bh; nt001 = nt001./bh; nt101 = nt101./bh; 
nt011 = nt011./bh;  nt103 = nt103./bh; nt013 = nt013./bh; 
nt203 = nt203./bh; nt023 = nt023./bh; nt113 = nt113./bh;
nt303 = nt303./bh; nt033 = nt033./bh; 
nt213 = nt213./bh; nt123 = nt123./bh;


%errors (comparison to numerical integration)
%these are really errors from numerical integration, not from our
%computation
errorsT003 = abs(t003 - nt003);
errorsT001 = abs(t001 - nt001);
errorsT101 = abs(t101 - nt101);
errorsT011 = abs(t011 - nt011);
errorsT103 = abs(t103 - nt103);
errorsT013 = abs(t013 - nt013); 
errorsT203 = abs(t203 - nt203);
errorsT023 = abs(t023 - nt023);
errorsT113 = abs(t113 - nt113);
errorsT303 = abs(t303 - nt303);
errorsT033 = abs(t033 - nt033);
errorsT213 = abs(t213 - nt213); 
errorsT123 = abs(t123 - nt123);
semilogy(errorsT003,'-o'); hold on; semilogy(errorsT001,'-o');
semilogy(errorsT101); semilogy(errorsT011,'-o'); 
semilogy(errorsT103, 'o-'); semilogy(errorsT203, 'o-'); 
semilogy(errorsT023, 'o-'); semilogy(errorsT113, 'o-');
semilogy(errorsT303, 'o-'); semilogy(errorsT033, 'o-'); 
semilogy(errorsT213, 'o-'); semilogy(errorsT123, 'o-'); 
hold off

function out=nIntegrateLineAlphaRp1(x0DotV,ell1,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a)(x0DotV(i)+a.*ell1).^2 + R0(i).^2 - (x0DotV(i)).^2;
    Rp1=@(a) a.* R2(a).^(1/2);
    out(i) = ell1.*integral(Rp1,0,1);
end
end

function out=nIntegrateLineRm1(x0DotV,ell1,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a)(x0DotV(i)+a.*ell1).^2 + R0(i).^2 - (x0DotV(i)).^2;
    Rm1=@(a) R2(a).^(-1/2);
    out(i) = ell1.*integral(Rm1,0,1);
end
end

function out=nIntegrateRm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    Rm3=@(a,b) R2(a,b).^(-3/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(Rm3,0,1,0,bmax);
end
end

function out=nIntegrateRm1(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    Rm1=@(a,b) R2(a,b).^(-1/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(Rm1,0,1,0,bmax);
end
end

function out=nIntegrateAlphaRm1(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    alphaRm1=@(a,b) a .* R2(a,b).^(-1/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(alphaRm1,0,1,0,bmax);
end%compare with numerical integration
end

function out=nIntegrateBetaRm1(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    betaRm1=@(a,b) b .* R2(a,b).^(-1/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(betaRm1,0,1,0,bmax);
end
end

function out=nIntegrateAlphaRm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    alphaRm3=@(a,b) a .* R2(a,b).^(-3/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(alphaRm3,0,1,0,bmax);
end%compare with numerical integration
end

function out=nIntegrateBetaRm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    alphaRm3=@(a,b) b .* R2(a,b).^(-3/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(alphaRm3,0,1,0,bmax);
end%compare with numerical integration
end

function out=nIntegrateAlpha2Rm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    alphaRm3=@(a,b) a.^2 .* R2(a,b).^(-3/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(alphaRm3,0,1,0,bmax);
end%compare with numerical integration
end

function out=nIntegrateBeta2Rm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    alphaRm3=@(a,b) b.^2 .* R2(a,b).^(-3/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(alphaRm3,0,1,0,bmax);
end%compare with numerical integration
end

function out=nIntegrateAlphaBetaRm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    alphaRm3=@(a,b) a.*b .* R2(a,b).^(-3/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(alphaRm3,0,1,0,bmax);
end%compare with numerical integration
end

function out=nIntegrateAlpha3Rm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    alpha3Rm3=@(a,b) a.^3 .* R2(a,b).^(-3/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(alpha3Rm3,0,1,0,bmax);
end
end

function out=nIntegrateBeta3Rm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    beta3Rm3=@(a,b) b.^3 .* R2(a,b).^(-3/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(beta3Rm3,0,1,0,bmax);
end
end

function out=nIntegrateAlpha2BetaRm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    beta3Rm3=@(a,b) a.^2 .* b .* R2(a,b).^(-3/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(beta3Rm3,0,1,0,bmax);
end
end

function out=nIntegrateAlphaBeta2Rm3(x0DotV,x0DotW,vDotW,ell1,ell2,R0) 
%matlab doesn't support numerical integration on arrays so need for loop

loopMax = size(R0,1);
out = zeros(size(R0));

for i = 1:loopMax
    R2=@(a,b)(x0DotV(i)+a.*ell1).^2+(x0DotW(i)+b.*ell2).^2 ...
        +2.*ell1.*ell2.*a.*b.*vDotW + R0(i).^2-(x0DotV(i)).^2-(x0DotW(i)).^2;
    beta3Rm3=@(a,b) a .* b.^2 .* R2(a,b).^(-3/2);
    bmax=@(a)a;
    jacob = sqrt(1-vDotW.^2);
    out(i) = (ell1.*ell2.*jacob).*integral2(beta3Rm3,0,1,0,bmax);
end
end
