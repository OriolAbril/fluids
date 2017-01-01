% Given Ra_critic=1708 & K_critic=3.12 obtained in prac5a
% Let's plot the vertical velocity (its vertical profile) & Temperature

close all 
clear all
global Raoldns
Raoldns=30000;
set(0, 'DefaultAxesFontSize', 13)
% We supose sigma=0 and ky=0 as we have rolls (v=0)
d=5e-3; %m
K=1.4e-7; %diffusion constat
Gamma=10000.0; %temperature gradient
k=fminbnd(@minQ3a_ns,1,10); %critic wavenumber
Ra=minQ3a_ns(k); %Critic Rayleigh number
z=linspace(-2*pi/k,2*pi/k,250);
C=K/(Gamma*d^2);
lam=(Ra/k^4)^(1/3);
q0=k*(lam-1)^.5;
q=k*(1+lam/2*(1+1i*sqrt(3)))^.5;
Q=[1i*q0 -1i*q0 q -q conj(q) -conj(q)];
Q=Q([1 3 5]);
%M=[cosh(Q*z); Q.*sinh(Q*z); Q.^4.*cosh(Q*z)-2*k^2*Q.^2.*cosh(Q*z)+k^4*cosh(Q*z)];
%M=[cosh(Q/2); Q.*sinh(Q/2); Q.^4.*cosh(Q/2)-2*k^2*Q.^2.*cosh(Q/2)+k^4*cosh(Q/2)];
% We know that M(1,:)=W; M(2,:)=DW; M(3,:)=(D^2-k^2)^2 and
% W=A*M(1,1)+B*M(1,2)+C*M(1,3) so we hace to find the coefficients A,B,C
% As the system is over ranked we will impose A=1. To obtain B and C then
% we proceed in this way:

W=zeros(1,length(z)); DW=zeros(1,length(z));DDW=zeros(1,length(z));
ii=1;
for zz=z
M=[sinh(Q*zz); Q.*cosh(Q*zz); Q.^4.*sinh(Q*zz)-2*k^2*Q.^2.*sinh(Q*zz)+k^4*sinh(Q*zz)];
Mcoef=M(1:3,2:3);
b=[-M(1,1) -M(2,1) -M(3,1)]';
coef=Mcoef\b;
coef=[1; coef];
Wvec=sinh(Q(1)*zz)+coef(2)*sinh(Q(2)*zz)+coef(3)*sinh(Q(3)*zz);
W(1,ii)=Wvec;
DWvec=Q(1)*cosh(Q(1)*zz)+ coef(2)*Q(2)*cosh(Q(2)*zz)+coef(3)*Q(3)*cosh(Q(3)*zz);
DW(1,ii)=DWvec;
DDWvec=(Q(1)^4.*sinh(Q(1)*zz)-2*k^2*Q(1)^2.*sinh(Q(1)*zz)+k^4*sinh(Q(1)*zz))+coef(2)*(Q(2)^4.*sinh(Q(2)*zz)-2*k^2*Q(2)^2.*sinh(Q(2)*zz)+k^4*sinh(Q(2)*zz))...
    +coef(3)*(Q(3)^4.*sinh(Q(3)*zz)-2*k^2*Q(3)^2.*sinh(Q(3)*zz)+k^4*sinh(Q(3)*zz));
DDW(1,ii)=DDWvec;
ii=ii+1;
end
wpertvertaxis=C*real(W);
tpert=(1/(Ra*k^2))*real(DDW);
norm=max(tpert);
wpertaxisnorm=wpertvertaxis/norm;
tpertnorm=tpert/norm;
figure(1)
plot(z,wpertaxisnorm,'linewidth',2)
title('w perturbation')
xlabel('z')
ylabel('w')
grid on
figure(2)
plot(z,tpertnorm,'linewidth',2)
title('T perturbation')
xlabel('z')
ylabel('T')
grid on