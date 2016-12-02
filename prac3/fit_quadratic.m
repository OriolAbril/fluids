clear all
close all
Revec=4500:50:4950;
n=length(Revec);
Recrit=5772.22;
Revec=[Revec Recrit];
amp=zeros(1,n+1);
L0=6.1570;
L0name=num2str(L0+1e-12,'%5.2f');
ii=find(L0name == '.') ; L0name(ii)='p';
for k=1:n
    fileampli = ['SerTemp_L' L0name 'Re'  int2str(Revec(k))];
    load(fileampli)
    amp(k)=mean(avec(end-100:end));
end
figure(1)
plot(amp,Revec,'bo')
hold on
amesh=0:0.001:max(amp);
p1=polyfit(amp,Revec,4)
plot(amesh,polyval(p1,amesh),'r-')
figure(2)
plot(amp,Revec,'bo')
hold on
amesh=amesh';
Revec=Revec';
a123 = [amp.^4; amp.^3; amp.^2]'\(Revec-Recrit);
p2=[a123;0; Recrit]
plot(amesh,polyval(p2,amesh),'k-')
