function Rec= Re1critic(k,Re2)
global x W DX PHI PSI eta %Re1_in
f = @(Re1)tclinstabfun(Re1,Re2,k,x,W,DX,PHI,PSI);
options=optimoptions('fsolve', 'TolFun',1e-12,'TolX',1e-12);
Rec=fsolve(f,abs(Re2)*(1-sign(Re2))+Re2/eta*(1+sign(Re2))+20,options);

end

