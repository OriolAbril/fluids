close all
clear all
MM = 50;
eta=0.5;
[x,w,W,DX,DX2,PHI,DPHI,D2PHI,PSI,DPSI,D2PSI] = setupspec2ord(MM);
r1 = eta/(1-eta) ; r2 = 1/(1-eta); beta = 2/(r2-r1) ;
k = 3.16 ; Re1 = 68.19 ; Re2 = -100; 
z=linspace(0,4*pi/k,500);
[EVEC,eval]=tclinstabeigs(Re1,Re2,k,eta,x,W,DX,PHI,PSI);
PHI1 = [PSI;0*PSI]; PHI2 = [0*PSI;PHI];
r = r1 + (x+1)/beta;
[R,Z]=meshgrid(r,z);
Ur=real(PHI1*EVEC(1:MM,end));
Ur=Ur(1:length(r))/max(abs(Ur)); 
%Ur=Ur(1:length(r));
[UR,~]=meshgrid(Ur,z);
psi=R.*UR.*cos(k*Z);
figure(1)
contour(R,Z,psi,20)
colorbar()
xlabel('r')
ylabel('z')
figure(2)
contour3(R,Z,psi,40)
colorbar()
xlabel('r')
ylabel('z')
figure(3)
plot(r,Ur)
xlabel('r')
ylabel('u_r')
grid on


