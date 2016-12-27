% codi de la practica 5, part analitica
close all
clear all
set(0, 'DefaultAxesFontSize', 13)
global k
kvec=linspace(1.5,10,100);
Ravecns=kvec*0;
Ravecsf=kvec*0;
j=1;
Raoldns=40000;
Raoldsf=30000;
kvec=kvec(end:-1:1);
for k=kvec
    Raoldns=fsolve(@funQ3a_ns,Raoldns);
    Ravecns(j)=Raoldns;
    Raoldsf=fsolve(@funQ3a_sf,Raoldsf);
    Ravecsf(j)=Raoldsf;
    j=j+1;
end
figure(1)
plot(kvec,Ravecns)
title('No-slip B.C.')
xlabel('k')
ylabel('Ra')
axis([0 10 0 40000])
grid on
figure(2)
plot(kvec,Ravecsf)
title('Stress-free B.C.')
xlabel('k')
ylabel('Ra')
axis([0 10 0 40000])
grid on
figure(3)
plot(kvec,Ravecns,kvec,Ravecsf)
legend('No-slip','Stress-free')
title('Critical Ra number')
xlabel('k')
ylabel('Ra')
axis([0 10 0 40000])
grid on