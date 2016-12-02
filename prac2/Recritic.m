function Rec= Recritic(k,Re0)
f = @(x)lmax(x,k);
options=optimoptions('fsolve', 'TolFun',1e-12,'TolX',1e-12);
Rec=fsolve(f,Re0,options);

end

