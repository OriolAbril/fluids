function [x,w,W,DX,D2X,PHI,DPHI,D2PHI,PSI,DPSI,D2PSI] = setupspec2ord(M)
% X - Legendre setting 
N = ceil(3*M/2) ; if mod(N,2)==0 ; N=N+1; end
[x,w,DX] = glxdw(N) ; W = diag(w) ; ID = eye(N) ;
D2X = DX*DX ; 
% Legendre polynomials generator
L0 = 1 + 0*x ; L1 = x ;  LM = zeros(N,M) ; LM(:,1) = L0 ; LM(:,2) = L1 ;
for k = 1:M-2
    icol = k + 1 ;
    LM(:,icol+1) = (2*k+1)*x.*LM(:,icol)/(k+1) - (k/(k+1))*LM(:,icol-1) ;
end
% PHI-homogeneous and PSI basis
PHI = zeros(N,M) ; PSI = zeros(N,M) ; 
g0 = (1-x.^2) ; h0 = (1-x.^2).^2 ;
for icol = 1:M
    PHI(:,icol) = g0.*LM(:,icol);
    PSI(:,icol) = h0.*LM(:,icol);
end
DPHI = DX*PHI ; D2PHI = D2X*PHI; 
DPSI = DX*PSI; D2PSI = D2X*PSI; 


function [x,w,D] = glxdw(N)

% Computation of Gauss Legendre points & weights
for N = 1:N;
    beta = .5./sqrt(1-(2*(1:N-1)).^(-2));
    T = diag(beta,1) + diag(beta,-1); [V,D] = eig(T);
    xxx = diag(D);                       %  <- Gauss nodes
    www = 2*V(1,:).^2;                   %  <- Gauss weights
    [junk2,index2]=sort(xxx);
    x=xxx(index2(:));
    w=www(index2(:));
end

index = (1:N)';
D = zeros(N,N); a = zeros(N,1);
for k = 1:N
    notk = find(index~=k);
    a(k) = prod(x(k)-x(notk));
end
for k = 1:N
    notk = find(index~=k);
    D(notk,k) = (a(notk)/a(k))./(x(notk)-x(k));
    D(k,k) = sum(1./(x(k)-x(notk)));
end
