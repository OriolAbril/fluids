function detM = funQ3a_sf(Ra,k)
%stress free boundary condition
lam=(Ra/k^4)^(1/3);
q0=k*(lam-1)^.5;
q=k*(1+lam/2*(1+1i*sqrt(3)))^.5;
Q=[1i*q0 -1i*q0 q -q conj(q) -conj(q)];
Q=Q([1 3 5]);
zz=0.5;
M=[sinh(Q*zz); Q.*sinh(Q*zz); Q.^4.*sinh(Q*zz)-2*k^2*Q.^2.*sinh(Q*zz)+k^4*sinh(Q*zz)];
detM=abs(det(M));

end

