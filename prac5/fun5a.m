function [M,Q] = fun5a(k,Ra)
%genera M i Q per les condicions de contorn no slip
%q0=(-k^2+(2*k^6+k^2*Ra)^(1/3))^.5;
%q=(k^2+.5*(1+1i*3^.5)*(2*k^6+k^2*Ra)^(1/3))^.5;
lam=(Ra/k^4)^(1/3);
q0=k*(lam-1)^.5;
q=k*(1+lam/2*(1+1i*sqrt(3)))^.5;
Q=[1i*q0 -1i*q0 q -q conj(q) -conj(q)];
Q=Q([1 3 5]);
%M=[cosh(Q/2); Q.*sinh(Q/2); (Q.^2-k^2).^2.*cosh(Q/2)];
M=[cosh(Q/2); Q.*sinh(Q/2); Q.^4.*cosh(Q/2)-2*k^2*Q.^2.*cosh(Q/2)+k^4*cosh(Q/2)];

end

