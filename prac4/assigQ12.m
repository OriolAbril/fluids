%codi per trobar la corva d'estabilitat respecte Re1 i Re2
close all
clear all
MM = 24;
global eta x W DX PHI PSI %Re1_in
etavec=[0.1 0.5 0.9];
for eta=etavec
    [x,w,W,DX,DX2,PHI,DPHI,D2PHI,PSI,DPSI,D2PSI] = setupspec2ord(MM);
    k = 3.16 ; Re1 = 68.19 ; Re2 = 0;
    Re2vec=linspace(-150,150,200);
    Re1vec=Re2vec*0;
    kcritvec=Re1vec;
    j=1;
    kvec=0.7:0.1:5;
    diffvec=kvec;
    for Re2=Re2vec
        fun = @(k)Re1critic(k,Re2);
        kcrit = fminbnd(fun,0.7,8);
        kcritvec(j)=kcrit;
        Re1=Re1critic(kcrit,Re2);
        Re1vec(j)=Re1;
        j=j+1;
    end
    Re2mesh=0:0.05:max(Re2vec);
    Re1mesh=Re2mesh/eta;
    figure(2)
    plot(Re2vec,Re1vec)
    hold on
    figure(1)
    plot(Re2vec,kcritvec)
    xlabel('Re_2')
    ylabel('k')
    hold on
end
figure(2)
legend('\eta = 0.1','\eta = 0.5','\eta = 0.9')
xlabel('Re_2')
ylabel('Re_1')
figure(1)
legend('\eta = 0.1','\eta = 0.5','\eta = 0.9')
