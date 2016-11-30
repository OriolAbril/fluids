function dRe= diffRe(k)
h=1e-4;
Rep=Recritic(k+h,6000);
Rem=Recritic(k-h,6000);
dRe=(Rep-Rem)/(2.0*h);

end

