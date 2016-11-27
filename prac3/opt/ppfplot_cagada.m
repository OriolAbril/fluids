function ppfplot_cagada(vortfile,ampfile)
clf
close all

%load VORT_L12_Re4500_T500.mat
%load ampL12Re4500.mat
load(vortfile)
load(ampfile)





figure(1)
V = [-3.5:.5:3.5];V2 = [-3.5:.5:3.5];
x2 = linspace(0,max(max(X)),200);
y2 = linspace(min(min(Y)),max(max(Y)),120);
[X2,Y2] = ndgrid(x2,y2); 
VORTT2 = griddata(X,Y,VORTT,X2,Y2,'cubic');
[c,h]=contourf(X2,Y2,VORTT2,V) ;
%set(h, 'linestyle', 'none'); hold on
xlabel('x');ylabel('y','rotation',0)
axis equal ; axis([0 max(x2) min(y2) max(y2)]);
colormap(cool) ; colorbar
title('Vorticity (basic flow + perturbation)','fontsize',16)
set(gca,'fontsize',16)

figure(2)
for ii = 1:Lmax+1 ; 
    ccc = .75*rand(1,3) ;
    semilogy(tvec2,amps(ii,2:end),'color',ccc,'linewidth',2) ; hold on ;  
    text(tvec2(end),amps(ii,end),[int2str(Lmax+1-ii)],'color',ccc,'fontsize',16)
end
axis([0 tvec(end) 1e-16 .5])
title('2-Norm of Fourier Streamwise Modes','fontsize',16)
xlabel('Time','fontsize',16)
set(gca,'fontsize',16)

figure(3)
plot(tvec,abs(avec(2:end)),'-r','linewidth',2)
axis([0 tvec(end) 0 2]); grid on    
title('2-Norm of the perturbation field (basic flow not included)','fontsize',16)
xlabel('Time','fontsize',16)
set(gca,'fontsize',16)


