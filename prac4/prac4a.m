%codi per trobar la corva d'estabilitat respecte Re1 i Re2
close all
clear all
MM = 24;
global eta x W DX PHI PSI %Re1_in
eta=0.5;
[x,w,W,DX,DX2,PHI,DPHI,D2PHI,PSI,DPSI,D2PSI] = setupspec2ord(MM);
k = 3.16 ; Re1 = 68.19 ;
Re2vec=linspace(-200,200,200);
Re1vec=Re2vec*0;
j=1;
kvec=0.7:0.1:5;
diffvec=kvec;
for Re2=Re2vec
    fun = @(k)Re1critic(k,Re2);
    kcrit = fminbnd(fun,0.7,8);
    Re1=Re1critic(kcrit,Re2);
    Re1vec(j)=Re1;
    j=j+1;
end
Re2mesh=0:0.05:max(Re2vec);
Re1mesh=Re2mesh/eta;
figure(2)
plot(Re2vec,Re1vec,'k-')
xlabel('Re_2')
ylabel('Re_1')
hold on
plot(Re2mesh,Re1mesh,'r--')
grid on