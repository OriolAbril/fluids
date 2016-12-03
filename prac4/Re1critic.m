function Rec= Re1critic(k,Re2)
global x W DX PHI PSI %Re1_in
f = @(Re1)tclinstabfun(Re1,Re2,k,x,W,DX,PHI,PSI);
options=optimoptions('fsolve', 'TolFun',1e-12,'TolX',1e-12);
Rec=fsolve(f,2*abs(Re2)+10,options);

end

