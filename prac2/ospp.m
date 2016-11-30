% Orr-Sommerfeld operators function (plane Poiseuille 2D)
% A. Meseguer (2015) - PHYSFLU 
% Input: k0 - wavelength & Re - Reynolds number
% Output: EFUN/eval/y - eigenfunctions/values/Legendre nodes
  function [EFUN,eval,y] = ospp(k0,Re)
  M = 80; N = 3*M/2 ;  [y,w,DY] = glxdw(N) ; W = diag(w) ; ID = eye(N) ;
  D2Y = DY*DY ; D4Y = D2Y*D2Y ; 
  BF = diag(1-y.^2) ; DBF = diag(-2*y) ; D2BF = diag(-2*(1+0*y)); % Poiseuille
  % BF = diag(y) ; DBF = 1+0*diag(y) ; D2BF = 0*DBF; %Couette
  
% Legendre polynomials generator
  L0 = 1 + 0*y ; L1 = y ;  LM = zeros(N,M) ; LM(:,1) = L0 ; LM(:,2) = L1 ;
  for k = 1:M-2
   icol = k + 1 ;
   LM(:,icol+1) = (2*k+1)*y.*LM(:,icol)/(k+1) - (k/(k+1))*LM(:,icol-1) ;
 end

% PHI-homogeneous basis
  PHI = zeros(N,M) ;  h0 = (1-y.^2).^2 ;
  for icol = 1:M
   PHI(:,icol) = h0.*LM(:,icol);
   DPHI(:,icol) = DY*PHI(:,icol) ; 
   D2PHI(:,icol) = D2Y*PHI(:,icol) ;
  end

% Operators
  k02 = k0^2; k04 = k0^4 ;
  A = PHI'*W*(D2Y-k02*ID)*PHI; 
  B = PHI'*W*((1/Re)*(D2Y-k02*ID)^2 + 1i*k0*D2BF - 1i*k0*BF*(D2Y-k02*ID))*PHI  ;
  
% Spectra  
  [V,D]=eig(B,A); d = diag(D);
  [~, ii]=sort(-real(d)); V = V(:,ii); D = D(:,ii); d = d(ii);
  EFUN = PHI*V ; eval = d ;