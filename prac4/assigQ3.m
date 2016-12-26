close all
clear all
MM = 50;
eta=0.5;
[x,w,W,DX,DX2,PHI,DPHI,D2PHI,PSI,DPSI,D2PSI] = setupspec2ord(MM);
r1 = eta/(1-eta) ; r2 = 1/(1-eta); beta = 2/(r2-r1) ;
k = 3.16 ; Re1 = 68.19 ; Re2 = -100;
[EVEC,eval]=tclinstabeigs(Re1,Re2,k,eta,x,W,DX,PHI,PSI);
PHI1 = [PSI;0*PSI]; PHI2 = [0*PSI;PHI];
r = r1 + (x+1)/beta;
Ur=real(PHI1*EVEC(1:MM,end));
Ur=Ur(1:length(r))/max(abs(Ur));
Uteta=real(PHI2*EVEC(1:MM,end));
Uteta=Uteta(length(r)+1:end)/max(abs(Uteta));
figure(1)
plot(r,Ur,'linewidth',2)
xlabel('r')
ylabel('u_r')
grid on
figure(2)
plot(r,Uteta)
xlabel('r')
ylabel('u_{\theta}')
grid on
