function ppftstep_opt_save(Ain,Scal,Tvec,L0in,Rein)
%inicialitza tots els vectors i només guarda cada 100 tots els valors
% Channel flow (BDF4 + mAB4 or IMEX4)
    Time = Tvec(end)-Tvec(1); k0 = 2*pi/L0in ; dt = 0.01 ; iter = round(Time/dt); savetime=100;
    L = 21 ; M = 40 ; Lmax = L ; Lmax2 = (3*Lmax+1)/2 ; 
 
% Initial Condition
% Option I (random with scaled initial amplitude factor Scal)
   if Ain == 0
    a0 = zeros(2*Lmax+1,M); 
    for l = 0:Lmax
      pert = .5*ones(1,M).*10.^(-[0:M-1]-l);
      lM = l + Lmax + 1 ; a0(lM,:) = pert*Scal;
      lMc = -l + Lmax + 1 ; a0(lMc,:) = conj(a0(lM,:));
    end	    
   else
    load(Ain) ; a0 = a0*Scal;
   end	
   Re = Rein ; L0 = L0in ; k0 = 2*pi/L0 ; 

% Preparation: BD4 
  [DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,RHS,INVA,x,y] = setupBD4(Re,dt,L,M,k0);
% Preparation: RK4 
  [INVARK4,INVABRK4] = prepRK4(Re,L,M,k0,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,y);
  a1 = RK4step(dt,L,M,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,INVARK4,INVABRK4,a0);
  a2 = RK4step(dt,L,M,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,INVARK4,INVABRK4,a1);
  a3 = RK4step(dt,L,M,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,INVARK4,INVABRK4,a2);
% Nonlinear terms for mAB4:
  b0 = nonlin(a0,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,Lmax,Lmax2);
  b1 = nonlin(a1,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,Lmax,Lmax2);
  b2 = nonlin(a2,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,Lmax,Lmax2);
  b3 = nonlin(a3,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,Lmax,Lmax2);

  
    
% ***************************   TIME STEPPING ************************
% Time specs and initial condition
  t = 3*dt+Tvec(1) ; 
  
  [X,Y] = ndgrid(x,y); 

% Monitoring variables
  tvec=(0:3)*dt; 
  tvec2 = [tvec(1:4) (1:(round(iter/savetime)-1))*dt*savetime]+Tvec(1);
  tvec=tvec2;
  avec=zeros(1,round(iter/savetime)+3);
  avec(1:4) = [norm(a0) norm(a1) norm(a2) norm(a3)] ; 
  l = 5 ; lM =  l + Lmax + 1 ; 
  covec=zeros(1,round(iter/savetime)+3);
  covec(1:4) = [a0(lM,2) a1(lM,2) a2(lM,2) a3(lM,2)]; 
  l = 5 ; lM =  l + Lmax + 1 ; 
  coveceven=zeros(1,round(iter/savetime)+3);
  covecodd=zeros(1,round(iter/savetime)+3);
  coveceven(1:4) = [a0(lM,1) a1(lM,1) a2(lM,1) a3(lM,1)];
  covecodd(1:4) = [a0(lM,2) a1(lM,2) a2(lM,2) a3(lM,2)];

  amps = zeros(Lmax+1,round(iter/savetime)+3); ampaux = zeros(Lmax+1,1) ;
  for ii = 1:Lmax + 1; ampaux(ii) = norm(a0(ii,:)); end ; amps(:,1) = ampaux;
  for ii = 1:Lmax + 1; ampaux(ii) = norm(a1(ii,:)); end ; amps(:,2) = ampaux;
  for ii = 1:Lmax + 1; ampaux(ii) = norm(a2(ii,:)); end ; amps(:,3) = ampaux;
  for ii = 1:Lmax + 1; ampaux(ii) = norm(a3(ii,:)); end ; amps(:,4) = ampaux;

% Name of files where data will be written
  L0name=num2str(L0+1e-12,'%5.2f');
  ii=find(L0name == '.') ; L0name(ii)='p';
  fileampli = ['SerTemp_L' L0name 'Re'  int2str(Re)];  
  
% BD4mAB4  
  h = waitbar(0,'Progress');
  for itime = 3:iter
      waitbar(itime/iter);

   for l = 0:Lmax
    lM = l + Lmax + 1 ; 
    aaux0 = a0(lM,:).' ; baux0 = b0(lM,:).'; aaux1 = a1(lM,:).' ; baux1 = b1(lM,:).';
    aaux2 = a2(lM,:).' ; baux2 = b2(lM,:).'; aaux3 = a3(lM,:).' ; baux3 = b3(lM,:).';
    aaux1 = RHS(:,:,lM)*(48*aaux3-36*aaux2+16*aaux1- 3*aaux0 +...
                         dt*INVA(:,:,lM)*(48*baux3-72*baux2+48*baux1-12*baux0)) ;
    a4(lM,:) = aaux1.'; 
    lMc = -l + Lmax + 1 ;  a4(lMc,:) = conj(a4(lM,:));
   end

   a0 = a1 ; a1 = a2 ; a2 = a3 ; a3 = a4 ; 
   b0 = b1 ; b1 = b2 ; b2 = b3 ; t = t + dt ;

%  ************** Updating Nonlinear term  ***************************
   b3 = nonlin(a3,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,Lmax,Lmax2);
%  ***********************************************************

  
   if mod(itime,savetime) == 0  || itime == iter
    for ii = 1:Lmax + 1; ampaux(ii) = norm(a3(ii,:)); end ; 
      amps(:,itime/savetime+4) =  ampaux; 
      avec(itime/savetime+4) = norm(a3); 
      l = 1 ; lM =  l + Lmax + 1 ; 
      covec(itime/savetime+4) = a3(lM,2); 
      coveceven(itime/savetime+4) = a3(lM,1); 
      covecodd(itime/savetime+4) = a3(lM,2);
   end
   
   % Writing data
   % if mod(itime,5000) == 0 | itime == iter
   [ii1,ii2,foo]=find(abs(Tvec-t)<1e-6);
   if length(ii2) > 0
       %[t ii2]
       %pause
        save(fileampli,'tvec', 'tvec2', 'avec', 'covec',...
             'amps','Lmax','coveceven','covecodd')
        %end
   
   %fileampli = ['ampL' L0name 'Re'  int2str(Re)];
   
   %if mod(itime,5000) == 0 | itime == iter
       filestate = ['A_L' L0name '_Re' int2str(Re) '_T' int2str(round(t))];
       save(filestate,'a0','a1','a2','a3','Re','L','M','k0','dt','Lmax') 
       PSI =  real(EF2X*(PHI*a2.').'); PSIT = Y-(1/3)*Y.^3 + PSI;
       VORTT = D2X*PSIT + (D2Y*PSIT.').';
       filevort =  ['VORT_L' L0name '_Re' int2str(Re) '_T' int2str(round(t))];
       save(filevort,'X','Y', 'VORTT','Re','t','L0')
   end
  end
  close(h)
  return
  
% ***************************   SETUP RK4    *************************
   function [INVARK4,INVABRK4] = prepRK4(Re,L,M,k0,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,y)
  Lmax = L ; Lmax2 = (3*Lmax+1)/2 ;  N = 3*M/2 ; ID = eye(N) ;
  BF = diag(1-y.^2) ; DBF = diag(-2*y) ; D2BF = diag(-2*(1+0*y)) ;
  % time-stepping operators (time step set here)
  INVARK4 = zeros(M,M,Lmax+1) ; INVABRK4 = INVARK4 ; lcount = 0 ; ID2 = eye(M) ;
   for l=0:Lmax 
    lcount = lcount + 1 ; k = k0*l ; k2 = k^2; k4 = k^4 ;
    A = PHI'*W*(D2Y-k2*ID)*PHI; INVARK4(:,:,lcount) = inv(A); 
    B = PHI'*W*((1/Re)*(D2Y-k2*ID)^2 + i*k*D2BF - i*k*BF*(D2Y-k2*ID))*PHI ;
    INVABRK4(:,:,lcount) = inv(A)*B ; 
   end
   return
% ***************************   END SETUP RK4  *************************

% ***************************   RK4 STEP    *****************************************
 function ark1 = RK4step(dtrk4,L,M,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,INVARK4,INVABRK4,a0)
 Lmax = L ; Lmax2 = (3*Lmax+1)/2 ;  

% Time advancing (four stages)
  a1 = zeros(2*Lmax+1,M) ; a2 = a1 ; a3 = a1 ; 
  a = zeros(2*Lmax+1,M)  ; b = a ; c = a ; d = a ; 
  ark0 = a0 ;
  
  % 1st stage RK4
  b0 = nonlin(a0,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,Lmax,Lmax2);lcount = 0 ;
  for l = 0:Lmax
      lcount = lcount + 1 ;
      lM = l + Lmax + 1 ; aaux0 = a0(lM,:).';baux0=b0(lM,:).';
      aaux1 = dtrk4*(INVABRK4(:,:,lcount)*aaux0 + INVARK4(:,:,lcount)*baux0) ;
      a(lM,:) = aaux1.'; 
      lMc = -l + Lmax + 1 ;  a(lMc,:) = conj(a(lM,:));
  end
      
  % 2nd stage RK4
  b0 = nonlin(a0+0.5*a,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,Lmax,Lmax2);lcount = 0 ;
  for l = 0:Lmax
      lcount = lcount + 1 ;
      lM = l + Lmax + 1 ; aaux0 = a0(lM,:).'+ 0.5*a(lM,:).'; baux0 = b0(lM,:).';
      aaux1 = dtrk4*(INVABRK4(:,:,lcount)*aaux0 + INVARK4(:,:,lcount)*baux0) ;
      b(lM,:) = aaux1.'; 
      lMc = -l + Lmax + 1 ;  b(lMc,:) = conj(b(lM,:));
  end
  
  % 3rd stage RK4
  b0 = nonlin(a0+0.5*b,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,Lmax,Lmax2);lcount = 0 ;
  for l = 0:Lmax
      lcount = lcount + 1 ;
      lM = l + Lmax + 1 ; aaux0 = a0(lM,:).'+ 0.5*b(lM,:).'; baux0 = b0(lM,:).';
      aaux1 = dtrk4*(INVABRK4(:,:,lcount)*aaux0 + INVARK4(:,:,lcount)*baux0) ;
      c(lM,:) = aaux1.'; 
      lMc = -l + Lmax + 1 ;  c(lMc,:) = conj(c(lM,:));
  end
      
  % 4th stage RK4
  b0 = nonlin(a0+c,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,Lmax,Lmax2);lcount = 0 ;
  for l = 0:Lmax
      lcount = lcount + 1 ;
      lM = l + Lmax + 1 ; aaux0 = a0(lM,:).'+ c(lM,:).'; baux0 = b0(lM,:).';
      aaux1 = dtrk4*(INVABRK4(:,:,lcount)*aaux0 + INVARK4(:,:,lcount)*baux0) ;
      d(lM,:) = aaux1.'; 
      lMc = -l + Lmax + 1 ;  d(lMc,:) = conj(d(lM,:));
  end
 
  ark1 = a0 + (1/6)*(a+2*b+2*c+d); 
  return
% *************************** END  RK4 STEP ************************************

% *************************** SETUP BD4 ************************************
 function [DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,RHS,INVA,x,y] = setupBD4(Re,dt,L,M,k0)

% X - De-aliased Fourier setting  Lmax ODD
  G = 2*pi/k0 ; 
  Lmax = L ; Lmax2 = (3*Lmax+1)/2 ; NX = 3*Lmax+2 ; % Original
  x = ((0:NX-1)*G/NX)'; % Axial mesh
  EX2F = (1/NX)*exp(-1i*(2*pi/G)*(-Lmax2:Lmax2).'*x.') ; % From X to F (X2F)
  EF2X = exp(1i*(2*pi/G)*x*(-Lmax:Lmax)) ;               % From F to X (F2X)
  hx = 2*pi/NX ; column = [0 .5*(-1).^(1:NX-1)./sin((1:NX-1)*hx/2)];
  DX = (2*pi/G)*toeplitz(column,column([1 NX:-1:2])); D2X = DX^2 ;  

% Y - Legendre setting (exact quadratures if ....) 
  N = 3*M/2 ;  [y,w,DY] = glxdw(N) ; W = diag(w) ; ID = eye(N) ;
  D2Y = DY*DY ; D4Y = D2Y*D2Y ; 
  BF = diag(1-y.^2) ; DBF = diag(-2*y) ; D2BF = diag(-2*(1+0*y)) ;
  
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

% time-stepping operators (time step set here)
  RHS = zeros(M,M,2*Lmax+1) ; INVA = RHS ; lcount = 0 ; ID2 = eye(M) ;
   for l=-Lmax:Lmax % Complex conjugate not yet implemented (pending)
    lcount = lcount + 1 ; k = k0*l ; k2 = k^2; k4 = k^4 ;
    A = PHI'*W*(D2Y-k2*ID)*PHI; INVA(:,:,lcount) = inv(A); 
    B = PHI'*W*((1/Re)*(D2Y-k2*ID)^2 + 1i*k*D2BF - 1i*k*BF*(D2Y-k2*ID))*PHI ;
    LIN = inv(A)*B ; 
%    [VV,DD] = eig(LIN); dd = diag(DD) ; [foo,ii]=sort(real(dd));
%    dd = dd(ii) ; DD = DD(:,ii) ; VV = VV(:,ii) ;
%    if real(dd(end)) > 0; 
%     [l dd(end)]
%     pause 
%    end  
    
    RHS(:,:,lcount) = inv(25*ID2-dt*12*LIN);
   end
   return
% ***************************   END SETUP BD4   *************************

% ***************************      NONLIN  ******************************
function b = nonlin(a,DX,DY,D2X,D2Y,W,PHI,EF2X,EX2F,Lmax,Lmax2)
  F0 = real(EF2X*(PHI*a.').');
  LAPPSI = D2X*F0 + (D2Y*F0.').';
  NLT = (DX*F0).*(DY*LAPPSI.').' - ((DY*F0.').').*(DX*LAPPSI) ;
  %preb = EX2F*(NLT*Y2L'); b = preb(Lmax2+1-Lmax:Lmax2+1+Lmax,:) ;
  preb = EX2F*((NLT*W)*PHI); b = preb(Lmax2+1-Lmax:Lmax2+1+Lmax,:) ;
  return
% ***************************  END NONLIN  ******************************

% ***************************    GLXDW   ******************************
  function [x,w,D] = glxdw(N)
    for N = 1:N;
     beta = .5./sqrt(1-(2*(1:N-1)).^(-2));
     T = diag(beta,1) + diag(beta,-1); [V,D] = eig(T);
     xxx = diag(D);                       %  <- Gauss nodes
     www = 2*V(1,:).^2;                   %  <- Gauss weights
    [~,index2]=sort(xxx);
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
    return
% ***************************  END GLXDW   ******************************
    


