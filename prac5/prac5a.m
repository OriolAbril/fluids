% codi de la practica 5, part analitica
close all
clear all
global k
kvec=linspace(0.5,6,50);
Ravec=kvec*0;
j=1;
Raold=3000;
kvec=kvec(end:-1:1);
for k=kvec
    Raold=fsolve(@fun5c,Raold);
    Ravec(j)=Raold;
    [M,Q]=fun5a(k,Raold);
    j=j+1;
end
plot(kvec,Ravec)
axis([0 6 0 6000])
grid on
k=kvec;
Ra=2000;
lam=(Ra./k.^4).^(1/3);
figure(2)
plot(kvec,(-k.^2+(2*k.^6+k.^2*Ra).^(1/3)).^.5,k,k.*(lam-1).^.5)