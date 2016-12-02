function sigma=tclinstabfun(Re1,Re2,k,x,W,DX,PHI,PSI)
global eta
%eta = .5 ; %k = 3.16 ; Re1 = 68.19 ; Re2 = 0 ; 
r1 = eta/(1-eta) ; r2 = 1/(1-eta); beta = 2/(r2-r1) ; r = r1 + (x+1)/beta;
D = beta*DX ; DS = D + diag(1./r) ; 
A = (Re2-eta*Re1)/(1+eta); B = eta*(Re1-eta*Re2)/((1-eta)*(1-eta^2));
V = (A+B./r.^2) ;
L22 = D^2+diag(1./r)*D-diag(1./(beta*r.^2))-k^2*eye(length(x));L11 = L22^2 ;
L12 = -2*k^2*diag(V); L21 = -2*A*eye(length(x));
L = [L11 L12 ; L21  L22] ; 
M = [L22 zeros(length(x)) ; zeros(length(x)) eye(length(x))];
PHI1 = [PSI;0*PSI]; PHI2 = [0*PSI;PHI]; WW = [W 0*W ; 0*W W];

AA  = [PHI1'*(WW*(L*PHI1)) PHI1'*(WW*(L*PHI2)) ;...
       PHI2'*(WW*(L*PHI1)) PHI2'*(WW*(L*PHI2))];
BB  = [PHI1'*(WW*(M*PHI1)) PHI1'*(WW*(M*PHI2)) ;...
       PHI2'*(WW*(M*PHI1)) PHI2'*(WW*(M*PHI2))];

[~,EVAL] = eig(AA,BB); eval = diag(EVAL); sigma=max(real(eval));
return
