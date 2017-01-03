function kelvin_graphic(zeta,u,v,x,y,xp,yp,dxp,dyp,N,M,Np,Mp,titol)
%Kelvin_graphic funció per fer els plots del codi de kelvin i fer-lo més curt i
%entenedor
% Variables for plotting
zetap = zeros(N,M); zetap(:,:) = NaN;
up = zeros(Np,Mp); up(:,:) = NaN;
vp = zeros(Np,Mp); vp(:,:) = NaN;
% Plot zeta everywhere except at the ghost points
for i=2:M-1
    zetap(2:N-1,i)=zeta(i,2:N-1);
end
% Calculate uc vc at the cells center and only every dp points
for i=2:dxp:M-1
    ii=1+(i-2)/dxp;
    up(1:Np,ii)=(u(i,2:dyp:N-1)+u(i+1,2:dyp:N-1))/2;
end
for j=2:dyp:N-1
    jj=1+(j-2)/dyp;
    vp(jj,1:Mp)=(v(2:dxp:M-1,j)+v(2:dxp:M-1,j+1))/2;
end

% Show initial free surface and velocity in a contour plot
pcolor(x,y,zetap)
shading flat
colorbar
hold on
quiver(xp,yp,up,vp,0.8,'k');
hold off
xlabel('Alongshore coordinate (m)')
ylabel('Cross-shore coordinate (m)')
title(titol)




end

