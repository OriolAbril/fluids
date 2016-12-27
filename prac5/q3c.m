%% No slip boundary conditions
close all 
clear all
set(0, 'DefaultAxesFontSize', 13)

% We supose sigma=0 and ky=0 as we have rolls (v=0)
d=5e-3; % m
K=1.4e-7; % diffusion constat
Gamma=10000.0; % temperature gradient
k=3.12; % wavenumber
Ra=1708; % Rayleigh number
X=linspace(0,2,100);
Z=linspace(0,2*pi/k,100);
[x,z]=meshgrid(X,Z);
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

W=zeros(length(X)); DW=zeros(length(X));DDW=zeros(length(X));
ii=1;
for zz=0:4*pi/((length(Z)-1)*k):4*pi/k
M=[sinh(Q*zz); Q.*cosh(Q*zz); Q.^4.*sinh(Q*zz)-2*k^2*Q.^2.*sinh(Q*zz)+k^4*sinh(Q*zz)];
Mcoef=M(1:3,2:3);
b=[-M(1,1) -M(2,1) -M(3,1)]';
coef=Mcoef\b;
coef=[1; coef];
Wvec=sinh(Q(1)*zz)+coef(2)*sinh(Q(2)*zz)+coef(3)*sinh(Q(3)*zz);
W(ii,:)=ones(1,length(X))* Wvec;
DWvec=Q(1)*cosh(Q(1)*zz)+ coef(2)*Q(2)*cosh(Q(2)*zz)+coef(3)*Q(3)*cosh(Q(3)*zz);
DW(ii,:)=ones(1,length(X))* DWvec;
DDWvec=(Q(1)^4.*sinh(Q(1)*zz)-2*k^2*Q(1)^2.*sinh(Q(1)*zz)+k^4*sinh(Q(1)*zz))+coef(2)*(Q(2)^4.*sinh(Q(2)*zz)-2*k^2*Q(2)^2.*sinh(Q(2)*zz)+k^4*sinh(Q(2)*zz))...
    +coef(3)*(Q(3)^4.*sinh(Q(3)*zz)-2*k^2*Q(3)^2.*sinh(Q(3)*zz)+k^4*sinh(Q(3)*zz));
DDW(ii,:)=ones(1,length(X))* DDWvec;
ii=ii+1;
end


wpert=C*real(W.*(cos(k*x)+1i*sin(k*x)));

tpert=(1/(Ra*k^2))*real(DDW.*(cos(k*x)+1i*sin(k*x)));

upert=(C/k)*real(DW.*(1i*cos(k*x)-sin(k*x)));

figure(1)
contour3(x,z,wpert,100)
title('w perturbation no-slip bc')
xlabel('x')
ylabel('z')
zlabel('w')
figure(2)
contour3(x,z,tpert,100)
title('T perturbation no-slip bc')
xlabel('x')
ylabel('z')
zlabel('T')
figure(3)
contour(x,z,upert,70)
title('u perturbation no-slip bc')
xlabel('x')
ylabel('z')
hold on
quiver(x,z,upert,wpert)
figure(4)
contour(x,z,upert,70)
title('u perturbation no-slip bc')
xlabel('x')
ylabel('z')
hold on
for i = 1:numel(x)
    mody=sqrt(upert(i)^2+wpert(i)^2); 
    %normalitzem les fletxes de manera que totes tinguin modul 1
    upert(i) = upert(i)/mody;
    wpert(i) = wpert(i)/mody;
end
quiver(x,z,upert,wpert)

%% stress-free boundary conditions
W=zeros(length(X)); D2W=zeros(length(X));DDW=zeros(length(X));
ii=1;
for zz=0:4*pi/(length(Z)-1)*k):4*pi/k
M=[sinh(Q*zz); (Q.^2).*sinh(Q*zz); Q.^4.*sinh(Q*zz)-2*k^2*Q.^2.*sinh(Q*zz)+k^4*sinh(Q*zz)];
Mcoef=M(1:3,2:3);
b=[-M(1,1) -M(2,1) -M(3,1)]';
coef=Mcoef\b;
coef=[1; coef];
Wvec=sinh(Q(1)*zz)+coef(2)*sinh(Q(2)*zz)+coef(3)*sinh(Q(3)*zz);
W(ii,:)=ones(1,length(X))* Wvec;
D2Wvec=Q(1)*sinh(Q(1)*zz)+ coef(2)*Q(2)*sinh(Q(2)*zz)+coef(3)*Q(3)*sinh(Q(3)*zz);
D2W(ii,:)=ones(1,length(X))* D2Wvec;
DDWvec=(Q(1)^4.*sinh(Q(1)*zz)-2*k^2*Q(1)^2.*sinh(Q(1)*zz)+k^4*sinh(Q(1)*zz))+coef(2)*(Q(2)^4.*sinh(Q(2)*zz)-2*k^2*Q(2)^2.*sinh(Q(2)*zz)+k^4*sinh(Q(2)*zz))...
    +coef(3)*(Q(3)^4.*sinh(Q(3)*zz)-2*k^2*Q(3)^2.*sinh(Q(3)*zz)+k^4*sinh(Q(3)*zz));
DDW(ii,:)=ones(1,length(X))* DDWvec;
ii=ii+1;
end

wpert=C*real(W.*(cos(k*x)+1i*sin(k*x)));

tpert=(1/(Ra*k^2))*real(DDW.*(cos(k*x)+1i*sin(k*x)));

upert=(C/k)*real(D2W.*(1i*cos(k*x)-sin(k*x)));

figure(5)
contour3(x,z,wpert,100)
title('w perturbation stress-free bc')
xlabel('x')
ylabel('z')
zlabel('w')
figure(6)
contour3(x,z,tpert,100)
title('T perturbation stress-free bc')
xlabel('x')
ylabel('z')
zlabel('T')
figure(7)
contour(x,z,upert,70)
title('u perturbation stress-free bc')
xlabel('x')
ylabel('z')
hold on
quiver(x,z,upert,wpert)
figure(8)
contour(x,z,upert,70)
title('u perturbation stress-free bc')
xlabel('x')
ylabel('z')
hold on
for i = 1:numel(x)
    mody=sqrt(upert(i)^2+wpert(i)^2); 
    %normalitzem les fletxes de manera que totes tinguin modul 1
    upert(i) = upert(i)/mody;
    wpert(i) = wpert(i)/mody;
end
quiver(x,z,upert,wpert)