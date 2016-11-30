options=optimoptions('fsolve', 'TolFun',1e-12,'TolX',1e-12);
kk=linspace(0.9,1.07,20);
dRe=kk;
ii=1;
%for k=kk
%    dRe(ii)=diffRe(k);
%    ii=ii+1;
%end
plot(kk,dRe)
kcrit=secante(0.95,1.05,1e-12,2000,@diffRe)
Recritic(kcrit,6000)
