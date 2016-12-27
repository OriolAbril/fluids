% codi de la practica 5, part analitica
close all
clear all
set(0, 'DefaultAxesFontSize', 13)
global Raoldns Raoldsf
kvec=linspace(1.5,10,100);
Ravecns=kvec*0;
Ravecsf=kvec*0;
j=1;
Raoldns=40000;
Raoldsf=30000;
kvec=kvec(end:-1:1);
for k=kvec
    ns=@(Ra)funQ3a_ns(Ra,k);
    Raoldns=fsolve(ns,Raoldns);
    Ravecns(j)=Raoldns;
    sf=@(Ra)funQ3a_sf(Ra,k);
    Raoldsf=fsolve(sf,Raoldsf);
    Ravecsf(j)=Raoldsf;
    j=j+1;
end
kcritns=fminbnd(@minQ3a_ns,1,10);
kcritsf=fminbnd(@minQ3a_sf,1,10);
Racritns=minQ3a_ns(kcritns);
Racritsf=minQ3a_sf(kcritsf);
figure(1)
plot(kvec,Ravecns)
hold on
scatter(kcritns,Racritns)
title('No-slip B.C.')
xlabel('k')
ylabel('Ra')
axis([0 10 0 40000])
grid on
figure(2)
plot(kvec,Ravecsf)
hold on
scatter(kcritsf,Racritsf)
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